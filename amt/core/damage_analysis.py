import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union
import matplotlib.pyplot as plt
import pysam

# Import C++ extension
try:
    from amt.cpp_extensions import damage_patterns as dp_ext
except ImportError:
    print("Warning: C++ extensions not installed. Using slower Python implementations.")
    dp_ext = None

class DamageProfiler:
    """
    Analyzes damage patterns in ancient DNA sequencing data.
    
    This class identifies and quantifies characteristic damage patterns in ancient DNA,
    particularly focusing on C→T and G→A deamination patterns at the 5' and 3' ends
    of DNA fragments.
    """
    
    def __init__(self, min_length: int = 30, max_length: int = 100, 
                 window_size: int = 10, use_cpp_ext: bool = True):
        """
        Initialize the DamageProfiler.
        
        Args:
            min_length: Minimum fragment length to consider
            max_length: Maximum fragment length to consider
            window_size: Window size for damage pattern analysis at 5' and 3' ends
            use_cpp_ext: Whether to use C++ extensions for performance-critical operations
        """
        self.min_length = min_length
        self.max_length = max_length
        self.window_size = window_size
        self.use_cpp_ext = use_cpp_ext and dp_ext is not None
        self.damage_profile = None
        
    def analyze_bam(self, bam_path: Union[str, Path], 
                    reference_path: Union[str, Path],
                    threads: int = 1) -> Dict:
        """
        Analyze damage patterns from a BAM file aligned to a reference.
        
        Args:
            bam_path: Path to the BAM file
            reference_path: Path to the reference genome
            threads: Number of threads to use
            
        Returns:
            Dictionary containing damage statistics
        """
        bam_path = Path(bam_path)
        reference_path = Path(reference_path)
        
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_path}")
        if not reference_path.exists():
            raise FileNotFoundError(f"Reference file not found: {reference_path}")
        
        print(f"Analyzing damage patterns in {bam_path}...")
        
        # C++ extension
        if self.use_cpp_ext:
            damage_stats = dp_ext.analyze_damage(
                str(bam_path), 
                str(reference_path),
                self.min_length,
                self.max_length,
                self.window_size,
                threads
            )
            self.damage_profile = damage_stats
            return damage_stats
        
        return self._analyze_bam_python(bam_path, reference_path)
    
    def _analyze_bam_python(self, bam_path: Path, reference_path: Path) -> Dict:
        """Python implementation of damage pattern analysis"""
        damage_counts = {
            "5_prime": {base_from + base_to: np.zeros(self.window_size) 
                         for base_from in "ACGT" for base_to in "ACGT" if base_from != base_to},
            "3_prime": {base_from + base_to: np.zeros(self.window_size)
                         for base_from in "ACGT" for base_to in "ACGT" if base_from != base_to},
            "coverage": {
                "5_prime": np.zeros(self.window_size),
                "3_prime": np.zeros(self.window_size)
            },
            "fragment_length": np.zeros(self.max_length + 1)
        }
        
        with pysam.AlignmentFile(bam_path, "rb") as bam, \
             pysam.FastaFile(reference_path) as ref:
            
            for read in bam.fetch():
                if read.is_unmapped or read.is_secondary or read.is_duplicate:
                    continue
                
                if read.query_length < self.min_length or read.query_length > self.max_length:
                    continue
                
                # Update fragment length distribution
                damage_counts["fragment_length"][read.query_length] += 1
                
                # Get reference sequence for this alignment
                ref_seq = ref.fetch(read.reference_name, read.reference_start, read.reference_end)
                
                # Analyze mismatches at 5' and 3' ends
                for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    ref_base = ref_seq[ref_pos - read.reference_start].upper()
                    read_base = read.query_sequence[read_pos].upper()
                    
                    # Skip if bases match or if either is N
                    if ref_base == read_base or ref_base == 'N' or read_base == 'N':
                        continue
                    
                    # Check if we're in the 5' window
                    if read_pos < self.window_size:
                        damage_counts["5_prime"].get(ref_base + read_base, np.zeros(self.window_size))[read_pos] += 1
                        damage_counts["coverage"]["5_prime"][read_pos] += 1
                    
                    # Check if we're in the 3' window
                    if read_pos >= read.query_length - self.window_size:
                        pos_from_end = read.query_length - read_pos - 1
                        damage_counts["3_prime"].get(ref_base + read_base, np.zeros(self.window_size))[pos_from_end] += 1
                        damage_counts["coverage"]["3_prime"][pos_from_end] += 1
        
        # Calculate damage frequencies
        damage_freq = {
            "5_prime": {},
            "3_prime": {},
            "fragment_length_dist": damage_counts["fragment_length"] / np.sum(damage_counts["fragment_length"])
        }
        
        for end in ["5_prime", "3_prime"]:
            for pattern, counts in damage_counts[end].items():
                if pattern in ["C→T", "G→A"]:  # Focus on characteristic aDNA damage
                    coverage = damage_counts["coverage"][end]
                    # Avoid division by zero
                    coverage_mask = coverage > 0
                    freq = np.zeros_like(counts)
                    freq[coverage_mask] = counts[coverage_mask] / coverage[coverage_mask]
                    damage_freq[end][pattern] = freq
        
        self.damage_profile = {
            "counts": damage_counts,
            "frequencies": damage_freq
        }
        
        return self.damage_profile
    
    def plot_damage_patterns(self, output_dir: Union[str, Path] = None) -> Dict[str, plt.Figure]:
        """
        Generate plots visualizing the damage patterns.
        
        Args:
            output_dir: Directory to save plots to (optional)
            
        Returns:
            Dictionary of matplotlib figures
        """
        if self.damage_profile is None:
            raise ValueError("No damage profile available. Run analyze_bam first.")
            
        output_dir = Path(output_dir) if output_dir else None
        if output_dir:
            output_dir.mkdir(parents=True, exist_ok=True)
            
        figures = {}
        
        # Plot C→T and G→A misincorporation frequencies
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        
        # 5' end
        x = np.arange(1, self.window_size + 1)
        if "C→T" in self.damage_profile["frequencies"]["5_prime"]:
            ax[0].plot(x, self.damage_profile["frequencies"]["5_prime"]["C→T"], 
                      'r-o', label='C→T')
        if "G→A" in self.damage_profile["frequencies"]["5_prime"]:
            ax[0].plot(x, self.damage_profile["frequencies"]["5_prime"]["G→A"], 
                      'b-o', label='G→A')
        ax[0].set_title("5' Damage Pattern")
        ax[0].set_xlabel("Position from 5' end")
        ax[0].set_ylabel("Frequency")
        ax[0].legend()
        ax[0].grid(True, alpha=0.3)
        
        # 3' end
        if "C→T" in self.damage_profile["frequencies"]["3_prime"]:
            ax[1].plot(x, self.damage_profile["frequencies"]["3_prime"]["C→T"], 
                      'r-o', label='C→T')
        if "G→A" in self.damage_profile["frequencies"]["3_prime"]:
            ax[1].plot(x, self.damage_profile["frequencies"]["3_prime"]["G→A"], 
                      'b-o', label='G→A')
        ax[1].set_title("3' Damage Pattern")
        ax[1].set_xlabel("Position from 3' end")
        ax[1].set_ylabel("Frequency")
        ax[1].legend()
        ax[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        figures["damage_patterns"] = fig
        
        if output_dir:
            fig.savefig(output_dir / "damage_patterns.png", dpi=300)
            
        # Plot fragment length distribution
        fig, ax = plt.subplots(figsize=(10, 5))
        x = np.arange(self.min_length, self.max_length + 1)
        ax.bar(x, self.damage_profile["frequencies"]["fragment_length_dist"][self.min_length:])
        ax.set_title("Fragment Length Distribution")
        ax.set_xlabel("Fragment Length (bp)")
        ax.set_ylabel("Frequency")
        ax.grid(True, alpha=0.3)
        
        figures["fragment_length"] = fig
        
        if output_dir:
            fig.savefig(output_dir / "fragment_length.png", dpi=300)
            
        return figures

def calculate_authentication_score(damage_profile: Dict, 
                                  min_ct_freq: float = 0.05) -> float:
    """
    Calculate an authentication score based on damage patterns.
    
    Args:
        damage_profile: Damage profile from DamageProfiler
        min_ct_freq: Minimum C→T frequency at 5' end to consider authentic
        
    Returns:
        Authentication score (0-1)
    """
    if not damage_profile:
        return 0.0
    
    # Check if C→T deamination pattern is present at 5' end
    ct_5prime = damage_profile["frequencies"]["5_prime"].get("C→T", None)
    if ct_5prime is None or len(ct_5prime) == 0:
        return 0.0
    
    # Calculate score based on C→T frequency at position 1
    max_ct_freq = ct_5prime[0]
    
    # Sigmoidal scaling between 0 and 1
    if max_ct_freq < min_ct_freq:
        return 0.0
    elif max_ct_freq > 0.25:  # Strong damage signal
        return 1.0
    else:
        return (max_ct_freq - min_ct_freq) / (0.25 - min_ct_freq)