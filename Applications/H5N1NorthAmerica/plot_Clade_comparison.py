import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import pandas as pd
import os
import random
import sys
import subprocess
from datetime import datetime, timedelta

sys.path.append('/Users/nmueller/Documents/github/CoInfection-Material/Applications/NetworkViz/')
import baltic_bacter as bt


def setup_matplotlib():
    """Set up matplotlib defaults"""
    typeface = 'helvetica'
    mpl.rcParams['font.weight'] = 300
    mpl.rcParams['axes.labelweight'] = 300
    mpl.rcParams['font.family'] = typeface
    mpl.rcParams['font.size'] = 12
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 12
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['legend.fontsize'] = 10


def define_constants():
    """Define all constants used in the analysis"""
    clades = ["B3.13", "D1.1"]
    rate_shift_str = '0 1.25 2.5 3.75 5 6 50 200'
    rate_shifts = np.array([float(x) for x in rate_shift_str.split()])
    
    # Set random seed for reproducibility
    random.seed(6465546)
    np.random.seed(6465546)
    
    # Define the MRSI (most recent sample isolation date)
    mrsi = datetime(2025, 2, 17)
    mrsi_hpai = mrsi
    mrsi_lpai = datetime(2024, 8, 22)
    
    # Define segment order
    segment_order = ["HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA"]
    
    # Color definitions
    methods_colors = {
        "skygrowth": "#E41A1C",      # red
        "skyline": "#377EB8",         # blue
        "skyline_Ne": "#4DAF4A"       # green
    }
    
    species_colors = {
        "cow": "#238B45",             # teal
        "bird": "#D95F02"              # burnt orange
    }
    
    clade_colors = {
        "B3.13": "#238B45",           # teal
        "D1.1": "#D95F02"              # burnt orange
    }
    
    lineage_colors = {
        "HPAI": "#d62728",            # red (matching plot_glm.py)
        "LPAI": "#2ca02c",            # green (matching plot_glm.py)
        "unknown": "#4DAF4A"          # green
    }
    
    path = '/Users/nmueller/Documents/github/CoInfection-Material/Applications/H5N1NorthAmerica/'
    
    return (clades, rate_shifts, mrsi, mrsi_hpai, mrsi_lpai, segment_order,
            methods_colors, species_colors, clade_colors, lineage_colors, path)


def run_beast_commands(path, force=False):
    """Run BEAST commands to process the data"""
    beast_path = "/Applications/BEAST\\ 2.7.7/bin/"
    combined_dir = os.path.join(path, 'combined')
    os.makedirs(combined_dir, exist_ok=True)
    
    if not force:
        # Check if main output files exist
        trees_output = os.path.join(path, 'combined/HLHxNx.skygrowth.trees')
        tree_output = os.path.join(path, 'combined/HLHxNx.skygrowth.tree')
        log_output = os.path.join(path, 'combined/HLHxNx.skygrowth.log')
        clades_trees = os.path.join(path, 'combined/HLHxNx.skygrowth.clades.trees')
        clades_tsv = os.path.join(path, 'combined/HLHxNx.skygrowth.clades.tsv')
        
        if all(os.path.exists(f) for f in [trees_output, tree_output, log_output, clades_trees, clades_tsv]):
            print("Main BEAST output files already exist. Skipping BEAST commands.")
            print("Use force=True to rerun all commands.")
            return
    
    # 1. Combine tree files
    print("Running BEAST logcombiner for trees...")
    logcombiner_cmd = f"{beast_path}logcombiner -burnin 20 -log ./out/HLHxNx.skygrowth.rep*.trees -o ./combined/HLHxNx.skygrowth.trees"
    subprocess.run(logcombiner_cmd, shell=True, cwd=path)
    
    # 2. Summarize network
    print("Running BEAST ReassortmentNetworkSummarize...")
    summarize_cmd = f"{beast_path}applauncher ReassortmentNetworkSummarize -burnin 0 -followSegment 0 -positions MCC ./combined/HLHxNx.skygrowth.trees ./combined/HLHxNx.skygrowth.tree"
    subprocess.run(summarize_cmd, shell=True, cwd=path)
    
    # 3. Combine log files
    print("Running BEAST logcombiner for logs...")
    logcombiner_log_cmd = f"{beast_path}logcombiner -burnin 20 -log ./out/HLHxNx.skygrowth.rep*.log -o ./combined/HLHxNx.skygrowth.log"
    subprocess.run(logcombiner_log_cmd, shell=True, cwd=path)
    
    # 4. Mark clades
    print("Running MarkCladesFromCladeFile...")
    mark_clades_cmd1 = f"{beast_path}applauncher MarkCladesFromCladeFile -burnin 0 -followSegment 0 -tree ./combined/HLHxNx.skygrowth.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.skygrowth.clades.trees"
    subprocess.run(mark_clades_cmd1, shell=True, cwd=path)
    
    mark_clades_cmd2 = f"{beast_path}applauncher MarkCladesFromCladeFile -burnin 0 -followSegment 0 -printTable true -tree ./combined/HLHxNx.skygrowth.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.skygrowth.clades.tsv"
    subprocess.run(mark_clades_cmd2, shell=True, cwd=path)
    
    # 5. Process segments
    segment_order = ["HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA"]
    for s, segment in enumerate(segment_order):
        print(f"Processing segment {segment}...")
        segment_cmd = f"{beast_path}applauncher MarkCladesFromCladeFile -burnin 0 -followSegment {s} -printSegment {s} -tree ./combined/HLHxNx.skygrowth.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HLHxNx.skygrowth.{segment}.trees"
        subprocess.run(segment_cmd, shell=True, cwd=path)
        
        treeannotator_cmd = f"{beast_path}treeannotator -burnin 0 -height keep ./combined/HLHxNx.skygrowth.{segment}.trees ./combined/HLHxNx.skygrowth.{segment}.tree"
        subprocess.run(treeannotator_cmd, shell=True, cwd=path)


def calculate_reassortment_rates(log_file, mrsi, path):
    """Calculate reassortment rate quantiles over time (same structure as plot_glm.py)"""
    # Read rate shifts from XML file (same as plot_glm.py)
    xml_file_path = os.path.join(path, 'xmls/HLHxNx.skygrowth.rep0.xml')
    if os.path.exists(xml_file_path):
        with open(xml_file_path, 'r') as xml_file:
            xml_content = xml_file.read()
        import re
        rate_shifts_match = re.findall(r'<stateNode id="rateShifts" spec="RealParameter" value="([^"]+)"/>', xml_content)
        if rate_shifts_match:
            rate_shifts = np.array([float(x) for x in rate_shifts_match[0].split()])
        else:
            # Fallback to hardcoded values if not found in XML
            rate_shift_str = '0 1.25 2.5 3.75 5 6 50 200'
            rate_shifts = np.array([float(x) for x in rate_shift_str.split()])
    else:
        # Fallback to hardcoded values if XML doesn't exist
        rate_shift_str = '0 1.25 2.5 3.75 5 6 50 200'
        rate_shifts = np.array([float(x) for x in rate_shift_str.split()])
    
    # Calculate time points (same as plot_glm.py)
    time_points = [mrsi - timedelta(days=shift * 365) for shift in rate_shifts]
    time_points_dec = [(t - datetime(1970, 1, 1)).days / 365.25 + 1970 for t in time_points]
    
    # Define quantiles for proper confidence intervals (same as plot_glm.py)
    quantiles_to_plot = [0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975]
    
    all_quantiles = []
    valid_times = []
    
    # Use 1-based indexing to match R code (i from 1 to len(rate_shifts))
    for i in range(1, len(rate_shifts) + 1):
        col_name = f"InfectedToRho.{i}"
        if col_name not in log_file.columns:
            continue
        
        rates_at_time = log_file[col_name].values
        
        # Calculate all quantiles for this time point (same as plot_glm.py)
        time_quantiles = []
        for q in quantiles_to_plot:
            time_quantiles.append(np.quantile(rates_at_time, q))
        all_quantiles.append(time_quantiles)
        
        # Use i-1 for rate_shifts array (0-based) but i for column name (1-based)
        valid_times.append(time_points_dec[i-1])
    
    # Convert to numpy arrays for easier handling (same as plot_glm.py)
    all_quantiles = np.array(all_quantiles)  # Shape: (time_points, quantiles)
    
    return valid_times, all_quantiles, quantiles_to_plot


def create_modified_tree_file(tree_path):
    """Create a modified copy of the tree file with HPAI+LPAI replaced with HPAI_LPAI for baltic parsing"""
    import re
    import tempfile
    
    # Create modified file path
    base_dir = os.path.dirname(tree_path)
    base_name = os.path.basename(tree_path)
    modified_path = os.path.join(base_dir, base_name.replace('.tree', '_modified.tree'))
    
    # Check if modified file already exists
    if os.path.exists(modified_path):
        return modified_path
    
    try:
        # Read original file
        with open(tree_path, 'r') as f:
            content = f.read()
        
        # Replace HPAI+LPAI with HPAI_LPAI so baltic_bacter can parse it
        modified_content = content.replace('HPAI+LPAI', 'HPAI_LPAI')
        
        # Write modified file
        with open(modified_path, 'w') as f:
            f.write(modified_content)
        
        return modified_path
    except Exception as e:
        print(f"Warning: Could not create modified tree file: {e}")
        return tree_path  # Fallback to original


def load_and_process_tree(path, segment="HA"):
    """Load and process the phylogenetic tree using baltic_bacter"""
    tree_file = f'combined/HLHxNx.skygrowth.{segment}.tree'
    tree_path = os.path.join(path, tree_file)
    
    if not os.path.exists(tree_path):
        print(f"Warning: Tree file not found: {tree_path}")
        return None
    
    # Create modified tree file with HPAI+LPAI -> HPAI_LPAI for HA segment only
    if segment == "HA":
        modified_path = create_modified_tree_file(tree_path)
    else:
        modified_path = tree_path
    
    print(f"Loading tree: {modified_path}")
    try:
        ll = bt.loadNexus(
            modified_path,
            date_fmt='%Y-%m-%d',
            treestring_regex=r'tree\s+[A-Za-z0-9_]+',
            verbose=False
        )
    except AssertionError:
        # Some tree files may include characters such as '.' or '-' in the name
        ll = bt.loadNexus(
            modified_path,
            date_fmt='%Y-%m-%d',
            treestring_regex=r'tree\s+[A-Za-z0-9_.-]+',
            verbose=False
        )
    ll.drawTree()
    
    return ll


def draw_ha_tree_with_clades(ax, ll, mrsi, mrsi_dec, lineage_colors, species_colors):
    """Draw the HA tree with clade coloring and cow/bird indicators"""
    fromval = 2020
    toval = mrsi_dec + 0.05
    
    # Add timeline shading (behind everything) - matching plot_glm.py
    for year in range(2020, 2026):
        ax.axvspan(year, year + 0.5, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)
    
    max_loc_seen = 0
    
    # Draw tree branches - time on x-axis, y-position on y-axis
    for k in ll.Objects:
        x = k.absoluteTime
        xp = k.parent.absoluteTime if k.parent else x
        y = k.y
        
        # Get posterior support for HPAI+LPAI reassortment from baltic traits
        # The modified tree file has HPAI_LPAI instead of HPAI+LPAI
        loc = 0
        if hasattr(k, 'traits'):
            # Try HPAI_LPAI first (from modified file)
            if 'HPAI_LPAI' in k.traits:
                try:
                    loc_val = k.traits['HPAI_LPAI']
                    if isinstance(loc_val, str):
                        loc = float(loc_val) if loc_val != 'NA' else 0
                    else:
                        loc = float(loc_val)
                except (ValueError, TypeError):
                    loc = 0
            # Fallback to other possible names
            elif 'HPAI+LPAI' in k.traits:
                try:
                    loc_val = k.traits['HPAI+LPAI']
                    if isinstance(loc_val, str):
                        loc = float(loc_val) if loc_val != 'NA' else 0
                    else:
                        loc = float(loc_val)
                except (ValueError, TypeError):
                    loc = 0
        
        # Clamp loc to [0, 1]
        if loc > 1:
            loc = 1
        if loc > max_loc_seen:
            max_loc_seen = loc
        
        # Color based on posterior support (gradient matching R code scale but with green)
        # R uses: scale_color_gradient2(low="black", mid="#377EB8", high="#377EB8", midpoint=0.7)
        # Using green (#2ca02c) instead of blue, same scale
        import matplotlib.colors as mcolors
        
        # Create gradient: black (0) -> green #2ca02c (0.7) -> green #2ca02c (1)
        # Green #2ca02c is RGB (44, 160, 44)
        if loc <= 0.7:
            # Interpolate from black (0,0,0) to green #2ca02c (44, 160, 44) for values 0-0.7
            t = loc / 0.7  # Normalize to 0-1
            r = int(t * 44)
            g = int(t * 160)
            b = int(t * 44)
            col = f'#{r:02x}{g:02x}{b:02x}'
        else:
            # Stay green for values > 0.7
            col = '#2ca02c'
        
        # Draw vertical branch (time direction)
        if not isinstance(k, bt.reticulation) and k.parent:
            ax.plot([x, xp], [y, y], color=col, linewidth=0.5, zorder=2,
                   solid_capstyle='round', solid_joinstyle='round')
        
        # Draw horizontal node connections
        if k.branchType == 'node' and len(k.children) >= 2:
            left = k.children[-1].y
            right = k.children[0].y
            ax.plot([x, x], [left, right], color=col, linewidth=0.5, zorder=2,
                   solid_capstyle='round', solid_joinstyle='round')
        
        # Draw tip points
        elif isinstance(k, bt.leaf):
            # Check if cow
            is_cow = False
            if hasattr(k, 'name'):
                is_cow = 'cow' in str(k.name).lower()
            
            # Outer circle (black)
            size_outer = 2.5 if is_cow else 1.0
            ax.scatter(x, y, s=size_outer*20, facecolor='black', edgecolor='none', zorder=4)
            
            # Inner circle (colored)
            size_inner = 1.5 if is_cow else 0.5
            # Use purple for cows to distinguish from green (branches/LPAI) and red (HPAI)
            tip_color = '#9467bd' if is_cow else '#999999'  # Purple for cows, gray for birds
            # Purple (cow) above grey (bird), but both below black outer circle
            tip_zorder = 5.5 if is_cow else 4.5
            ax.scatter(x, y, s=size_inner*20, facecolor=tip_color, edgecolor='none', zorder=tip_zorder)
    
    ax.set_xlim(fromval, toval)
    # Set y-axis limits to show entire tree (matching plot_glm.py approach)
    ax.set_ylim(-ll.ySpan * 0.05, ll.ySpan * 1.01)
    ax.set_xlabel('')  # Remove x-axis title
    ax.set_ylabel('')
    ax.set_xticks(range(2020, 2026))
    ax.set_xticklabels(range(2020, 2026))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    
    # Add legend for reassortment support and cow isolates (matching R version)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.cm as cm
    
    # Create legend elements for cow isolates
    legend_elements = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor='#999999', 
               markersize=6, markeredgecolor='black', markeredgewidth=0.5, label='bird'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor='#9467bd', 
               markersize=10, markeredgecolor='black', markeredgewidth=0.5, label='cow')
    ]
    
    # Position legend in bottom-left corner inside the figure
    # Use bbox_to_anchor with small offsets to ensure it's inside the plot area
    legend1 = ax.legend(handles=legend_elements, loc='lower left', fontsize=9, 
                        frameon=False, title='cow isolate', title_fontsize=9,
                        bbox_to_anchor=(0.02, 0.02))
    
    # Add colorbar for reassortment support (gradient matching R code scale but with green)
    # R uses: scale_color_gradient2(low="black", mid="#377EB8", high="#377EB8", midpoint=0.7)
    # Using green (#2ca02c) instead of blue, same scale
    colors_list = ['black', '#2ca02c']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('reassortment', colors_list, N=n_bins)
    
    # Create a scalar mappable for the colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    
    # Add colorbar inside the plot area - create axes manually for consistent PDF output
    # Position it higher (y=0.32) to avoid overlap with cow legend at y=0.02
    # Get the axes position in figure coordinates
    fig = ax.figure
    ax_pos = ax.get_position()
    
    # Calculate colorbar position in figure coordinates
    # Position in axes coordinates: x=0.02, y=0.32, width=0.15, height=0.03
    cbar_width = 0.15 * ax_pos.width
    cbar_height = 0.03 * ax_pos.height
    cbar_x = ax_pos.x0 + 0.02 * ax_pos.width
    cbar_y = ax_pos.y0 + 0.32 * ax_pos.height
    
    # Create colorbar axes in figure coordinates
    cax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
    cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label('posterior support for\nHPAI LPAI reassortment', 
                   fontsize=9, labelpad=5)
    cbar.set_ticks([0, 0.5, 1.0])
    cbar.set_ticklabels(['0', '0.5', '1+'])  # Matching R code breaks and labels
    cbar.ax.tick_params(labelsize=8)
    # Remove background from colorbar axes to ensure proper integration
    cbar.ax.patch.set_facecolor('none')
    cbar.ax.patch.set_edgecolor('none')
    # Ensure colorbar doesn't resize with figure
    cbar.ax.set_anchor('W')



def draw_segment_tree(ax, ll, mrsi_dec, lineage_colors, segment_name):
    """Draw a segment tree colored by lineage type"""
    fromval = 2020
    toval = mrsi_dec + 0.05
    
    # Add timeline shading (behind everything) - matching plot_glm.py
    for year in range(2020, 2026):
        ax.axvspan(year, year + 0.5, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)
    
    # Draw tree branches - time on x-axis, y-position on y-axis
    for k in ll.Objects:
        x = k.absoluteTime
        xp = k.parent.absoluteTime if k.parent else x
        y = k.y
        
        # Get lineage type
        lin_type = 'unknown'
        if hasattr(k, 'traits'):
            if 'type' in k.traits:
                lin_type = str(k.traits['type']).upper()
            # Also check for other possible trait names
            elif 'Lineage' in k.traits:
                lin_type = str(k.traits['Lineage']).upper()
        
        col = lineage_colors.get(lin_type, lineage_colors['unknown'])
        
        # Draw vertical branch (time direction)
        if not isinstance(k, bt.reticulation) and k.parent:
            ax.plot([x, xp], [y, y], color=col, linewidth=0.2, zorder=2,
                   solid_capstyle='round', solid_joinstyle='round')
        
        # Draw horizontal node connections
        if k.branchType == 'node' and len(k.children) >= 2:
            left = k.children[-1].y
            right = k.children[0].y
            ax.plot([x, x], [left, right], color=col, linewidth=0.2, zorder=2,
                   solid_capstyle='round', solid_joinstyle='round')
        
        # Draw tip points
        elif isinstance(k, bt.leaf):
            # Outer circle (black)
            ax.scatter(x, y, s=20, facecolor='black', edgecolor='none', zorder=4)
            # Inner circle (colored)
            ax.scatter(x, y, s=10, facecolor=col, edgecolor='none', zorder=5)
    
    ax.set_xlim(fromval, toval)
    # Set y-axis limits to show entire tree (matching plot_glm.py approach)
    ax.set_ylim(-ll.ySpan * 0.05, ll.ySpan * 1.01)
    ax.set_xlabel('')  # Remove x-axis title
    ax.set_ylabel('')
    ax.set_xticks(range(2020, 2026, 2))
    ax.set_xticklabels(range(2020, 2026, 2))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, labelleft=False)


def plot_reassortment_rates_skygrowth(ax, valid_times, all_quantiles, quantiles_to_plot, mrsi, fromval, toval):
    """Plot reassortment rates from skygrowth analysis (same structure as plot_glm.py)"""
    # Add timeline shading
    timewidth = 0.5
    for i in np.arange(fromval, toval, 2 * timewidth):
        ax.axvspan(i, i + timewidth, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)
    
    # Plot confidence intervals as shaded ribbons (same as plot_glm.py)
    confidence_intervals = [
        (0.95, 0.025, 0.975),  # 95% CI: 2.5% to 97.5%
        (0.90, 0.05, 0.95),    # 90% CI: 5% to 95%
        (0.80, 0.1, 0.9),      # 80% CI: 10% to 90%
        (0.60, 0.2, 0.8),      # 60% CI: 20% to 80%
        (0.40, 0.3, 0.7),      # 40% CI: 30% to 70%
        (0.20, 0.4, 0.6),      # 20% CI: 40% to 60%
    ]
    
    # Plot each confidence interval as a shaded ribbon (same as plot_glm.py)
    for ci_level, lower_q, upper_q in confidence_intervals:
        # Find the indices for these quantiles
        lower_idx = quantiles_to_plot.index(lower_q)
        upper_idx = quantiles_to_plot.index(upper_q)
        
        # Get the quantile values (values are already in log space)
        lower_values = all_quantiles[:, lower_idx]
        upper_values = all_quantiles[:, upper_idx]
        
        # Calculate alpha based on confidence level
        # Higher confidence = lighter shading, lower confidence = darker shading
        alpha = 0.1 + (1.0 - ci_level) * 0.4  # Range from 0.1 to 0.5
        
        # Create label for key confidence intervals
        label = f'{int(ci_level*100)}% CI' if ci_level in [0.95, 0.80, 0.40] else None
        
        # Create shaded ribbon for this confidence interval
        ax.fill_between(valid_times, lower_values, upper_values,
                       alpha=alpha, color='#C85A3C',  # Strong terracotta
                       label=label, zorder=2, linewidth=0)
    
    # Plot the median as a thick shaded ribbon (same as plot_glm.py)
    median_alpha = 0.6
    median_thickness = 0.01  # Small thickness for the median ribbon
    
    # Find median index (50% quantile)
    median_idx = quantiles_to_plot.index(0.5)
    median_values = all_quantiles[:, median_idx]
    
    # Create upper and lower bounds for median ribbon
    median_upper = median_values + median_thickness
    median_lower = median_values - median_thickness
    
    ax.fill_between(valid_times, median_lower, median_upper,
                   alpha=median_alpha, color='#C85A3C',  # Strong terracotta
                   label='Median', zorder=3, linewidth=0)
    
    ax.set_xlim(fromval, toval)
    ax.set_xlabel('')
    ax.set_ylabel('reassortment rate', labelpad=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Set x-axis ticks to only show full years
    year_ticks = np.arange(np.ceil(fromval), np.floor(toval) + 1, 1.0)
    ax.set_xticks(year_ticks)
    ax.set_xticklabels([int(year) for year in year_ticks])
    
    # Set logarithmic y-axis with custom breaks and labels (same as plot_glm.py)
    tick_positions = [np.log(0.025), np.log(0.05), np.log(0.1), np.log(0.2), np.log(0.4), 
                     np.log(0.8), np.log(1.6), np.log(3.2)]
    tick_labels = ["0.025", "0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "3.2"]
    
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)
    ax.set_ylim(np.log(0.0125), np.log(6.4))


def draw_network_tree(ax, ll, mrsi, mrsi_dec, colour_cycle, linewidth):
    """Draw the network tree showing reassortment events (similar to plot_glm.py)"""
    fromval = float(int(ll.root.absoluteTime + ll.treeHeight)) - 6
    toval = float(int(ll.root.absoluteTime + ll.treeHeight)) + 0.5
    timewidth = 0.5
    
    # Add timeline shading
    for i in np.arange(fromval, toval, 2 * timewidth):
        ax.axhspan(i, i + timewidth, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)
    
    # Initialize traits
    for k in ll.Objects:
        if 're' not in k.traits:
            k.traits['re'] = 0
    
    # Assign reassortment colors
    for k in sorted(ll.Objects, key=lambda w: w.height):
        if hasattr(k, 'contribution'):
            random_number = random.randint(0, len(colour_cycle) - 1)
            while random_number == k.traits['re']:
                random_number = random.randint(0, len(colour_cycle) - 1)
            
            if hasattr(k, 'children') and len(k.children) > 0:
                subtree = ll.traverse_tree(k.children[-1], include_condition=lambda w: True)
                for w in subtree:
                    w.traits['re'] = random_number
    
    # Draw tree branches (rotated 90 degrees like plot_glm.py)
    main_branch_color = '#2c3e50'
    secondary_branch_color = '#95a5a6'
    reassortment_color = '#e74c3c'
    leaf_outer_color = '#34495e'
    
    for k in ll.Objects:
        x = k.absoluteTime
        xp = k.parent.absoluteTime if k.parent else x
        if xp is not None:
            xp = max(xp, fromval + 0.000001)
        y = k.y
        col = colour_cycle[k.traits['re'] % len(colour_cycle)]
        
        if not isinstance(k, bt.reticulation):
            col_lin = secondary_branch_color
            lw_scale = 1.2
            if hasattr(k, 'traits') and k.traits.get('seg0') == 'true':
                col_lin = col
                lw_scale = 2.5
            
            # Rotated: horizontal branches become vertical
            ax.plot([y, y], [x, xp], color=col_lin, lw=linewidth * lw_scale,
                   solid_capstyle='round', solid_joinstyle='round', alpha=1, zorder=2)
        else:
            # Reassortment branches
            ax.plot([y, y], [x, xp], color=col, lw=linewidth * 1.5,
                   ls='--', solid_capstyle='round', solid_joinstyle='round',
                   alpha=1, zorder=1)
        
        if k.branchType == 'node' and len(k.children) >= 2:
            left = k.children[-1].y
            right = k.children[0].y
            
            col_lin1 = col
            lw_scale1 = 2.5
            col_lin2 = col
            lw_scale2 = 2.5
            
            # Rotated: vertical lines become horizontal
            ax.plot([left, k.y], [x, x], color=col_lin1, lw=linewidth * lw_scale1,
                   solid_capstyle='round', solid_joinstyle='round', alpha=1, zorder=2)
            ax.plot([k.y, right], [x, x], color=col_lin2, lw=linewidth * lw_scale2,
                   solid_capstyle='round', solid_joinstyle='round', alpha=1, zorder=2)
        
        elif isinstance(k, bt.leaf):
            # Rotated: leaf nodes
            ax.scatter(y, x, s=30, facecolor=leaf_outer_color, edgecolor='none', zorder=4)
            ax.scatter(y, x, s=15, facecolor=col, edgecolor='none', zorder=5)
        
        elif isinstance(k, bt.reticulation):
            # Reassortment nodes
            ax.scatter(k.target.y, x, s=20, facecolor=reassortment_color,
                      edgecolor='white', linewidth=0.5, zorder=4, alpha=0.9)
            ax.scatter(k.target.y, x, s=6, facecolor=reassortment_color, edgecolor='none', zorder=5)
            ax.plot([y, k.target.y], [x, x], color=col, lw=linewidth * 1.5,
                   ls='-', solid_capstyle='round', solid_joinstyle='round', alpha=1, zorder=1)
    
    # Finalize plot (rotated)
    ax.set_xticks([])
    ax.set_xlim(ll.ySpan * 1.01, -ll.ySpan * 0.05)  # Inverted x-axis
    ax.set_ylim(toval, fromval)  # toval (recent) at bottom, fromval (past) at top - inverted
    
    # Add y-axis with time labels
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_linewidth(0.5)
    ax.set_ylabel('Time (years)')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', labelsize=10)
    
    # Set y-axis ticks to show years
    year_ticks = np.arange(np.ceil(fromval), np.floor(toval) + 1, 1.0)
    ax.set_yticks(year_ticks)
    ax.set_yticklabels([int(year) for year in year_ticks])


def calculate_event_distribution(hpai_events, segment_order):
    """Calculate event distributions by lineage and event type"""
    lineage_type = ["HPAI", "LPAI"]
    event_types = ["HPAI", "HPAI+LPAI", "LPAI"]
    
    distr = []
    co_rea = []
    
    for sample in range(int(hpai_events['Sample'].min()), int(hpai_events['Sample'].max()) + 1):
        for l in lineage_type:
            n_events = hpai_events[(hpai_events['Sample'] == sample) & (hpai_events['Lineage'] == l)]
            
            for e in event_types:
                # Check if lineage is in event type (matches R's grepl logic)
                if l not in e:
                    continue
                
                n_events_e = n_events[n_events['Event'] == e]
                n_events_count = len(n_events_e)
                
                distr.append({
                    'Sample': sample,
                    'Lineage': l,
                    'Event': e,
                    'n_events': n_events_count
                })
                
                # Count segment reassortments (excluding HA, segment 0)
                # Loop over all segments not HA (starting from index 1)
                for j in range(1, len(segment_order)):
                    segment = segment_order[j]
                    segment_index = j  # Segments are encoded 1-7 for NA-PA
                    
                    count = 0
                    if 'Segments' in n_events_e.columns:
                        for seg_str in n_events_e['Segments'].values:
                            if pd.isna(seg_str):
                                continue
                            seg_str = str(seg_str).strip('{} ')
                            if not seg_str:
                                continue
                            try:
                                nums = [int(s.strip()) for s in seg_str.split(',')]
                            except ValueError:
                                continue
                            if segment_index in nums:
                                count += 1
                    
                    co_rea.append({
                        'Sample': sample,
                        'Lineage': l,
                        'Event': e,
                        'n_events': count,
                        'segment': segment
                    })
    
    return pd.DataFrame(distr), pd.DataFrame(co_rea)


def calculate_co_rea_quantiles(co_rea_df):
    """Calculate quantiles for co-reassortment counts"""
    co_rea_quantiles = []
    
    for s in co_rea_df['segment'].unique():
        for e in co_rea_df['Event'].unique():
            for l in co_rea_df['Lineage'].unique():
                counts = co_rea_df[(co_rea_df['segment'] == s) & 
                                  (co_rea_df['Event'] == e) & 
                                  (co_rea_df['Lineage'] == l)]['n_events'].values
                
                if len(counts) == 0:
                    continue
                
                q = 0.05
                lower = np.quantile(counts, q/2)
                upper = np.quantile(counts, 1 - q/2)
                mean_val = np.mean(counts)
                
                evname = "HPAI" if e == "HPAI" else "LPAI"
                
                co_rea_quantiles.append({
                    'segment': s,
                    'Event': e,
                    'Lineage': l,
                    'mean': mean_val,
                    'lower': lower,
                    'upper': upper,
                    'evname': evname
                })
    
    return pd.DataFrame(co_rea_quantiles)


def plot_event_distribution(ax, distr_df):
    """Plot violin plot of event distribution"""
    lineage_type = ["HPAI", "LPAI"]
    event_types = ["HPAI", "HPAI+LPAI", "LPAI"]
    
    # Create subplots for each lineage
    for idx, l in enumerate(lineage_type):
        data_subset = distr_df[distr_df['Lineage'] == l]
        
        if len(data_subset) == 0:
            continue
        
        # Filter events that contain the lineage
        valid_events = [e for e in event_types if l in e]
        data_subset = data_subset[data_subset['Event'].isin(valid_events)]
        
        # Create violin plot
        parts = ax.violinplot([data_subset[data_subset['Event'] == e]['n_events'].values 
                              for e in valid_events],
                             positions=range(len(valid_events)),
                             showmeans=False, showmedians=True)
        
        for pc in parts['bodies']:
            pc.set_facecolor('#4E79A7')
            pc.set_alpha(0.6)
        
        if 'cmedians' in parts:
            parts['cmedians'].set_color('black')
            parts['cmedians'].set_linewidth(1.5)
    
    ax.set_xlabel('Event type')
    ax.set_ylabel('Number of events')
    ax.set_xticks(range(len(event_types)))
    ax.set_xticklabels(event_types)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, axis='y')


def plot_co_reassortment(ax, co_rea_quantiles_df, lineage_colors, segment_order):
    """Plot co-reassortment events by segment"""
    # Filter for HPAI lineage only
    data = co_rea_quantiles_df[co_rea_quantiles_df['Lineage'] == 'HPAI'].copy()
    
    if len(data) == 0:
        return
    
    desired_order = segment_order[1:]  # skip HA
    segments = [seg for seg in desired_order if seg in set(data['segment'])]
    x_pos = np.arange(len(segments))
    width = 0.35
    dodge_offset = -0.3  # Match R's position_dodge(-0.3)
    
    for evname in ['HPAI', 'LPAI']:
        ev_data = data[data['evname'] == evname]
        if len(ev_data) == 0:
            continue
        
        means = []
        lowers = []
        uppers = []
        
        for s in segments:
            seg_data = ev_data[ev_data['segment'] == s]
            if len(seg_data) > 0:
                means.append(seg_data['mean'].values[0])
                lowers.append(seg_data['lower'].values[0])
                uppers.append(seg_data['upper'].values[0])
            else:
                means.append(0)
                lowers.append(0)
                uppers.append(0)
        
        offset = dodge_offset if evname == 'HPAI' else -dodge_offset
        color = lineage_colors[evname]
        
        # Plot points
        ax.scatter(x_pos + offset, means, color=color, s=50, zorder=3, 
                  label=evname)
        
        # Plot error bars
        yerr_lower = np.array(means) - np.array(lowers)
        yerr_upper = np.array(uppers) - np.array(means)
        ax.errorbar(x_pos + offset, means, yerr=[yerr_lower, yerr_upper],
                   fmt='none', color=color, capsize=3, capthick=1, 
                   linewidth=1.5, zorder=2)
    
    ax.set_xlabel('')  # Remove x-axis title
    ax.set_ylabel('events')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(segments)
    ax.legend(loc='upper left', fontsize=9, title='segment origin', frameon=False,
              bbox_to_anchor=(0, 1.18))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, axis='y')


def create_main_figure(path, mrsi, mrsi_dec, segment_order, lineage_colors, species_colors, 
                       hpai_events, distr_df, co_rea_df, co_rea_quantiles_df, 
                       valid_times, all_quantiles, quantiles_to_plot):
    """Create the main composite figure (Figure 2)"""
    fig = plt.figure(figsize=(12, 8), constrained_layout=True)
    # Layout: A (left, top), Reassortment rates (top right), B (left, bottom), C and D (bottom right, side by side)
    # Top row fills most of the height
    # Use constrained_layout to ensure all content fits within subplot boundaries
    gs = gridspec.GridSpec(3, 3, width_ratios=[3, 1, 1], height_ratios=[1, 3, 1], 
                          figure=fig)

    # Panel A: HA tree with clades (left, top 2 rows)
    ax_tree = fig.add_subplot(gs[0:2, 0])
    ll_ha = load_and_process_tree(path, "HA")
    if ll_ha:
        draw_ha_tree_with_clades(ax_tree, ll_ha, mrsi, mrsi_dec, lineage_colors, species_colors)
    ax_tree.text(-0.02, 1.02, 'A', transform=ax_tree.transAxes, 
                fontsize=16, fontweight='bold', va='bottom', ha='right')
    
    # Panel B: Reassortment rates (top right)
    ax_rates = fig.add_subplot(gs[0, 1:])
    if valid_times is not None and len(valid_times) > 0:
        fromval = mrsi_dec - 10  # Adjust as needed
        toval = mrsi_dec
        plot_reassortment_rates_skygrowth(ax_rates, valid_times, all_quantiles, quantiles_to_plot, 
                                         mrsi, fromval, toval)
    ax_rates.text(-0.02, 1.02, 'B', transform=ax_rates.transAxes, 
                   fontsize=16, fontweight='bold', va='bottom', ha='right')
    
    # Panel C: Co-reassortment plot (left, bottom row, below A)
    ax_corea = fig.add_subplot(gs[2, 0])
    plot_co_reassortment(ax_corea, co_rea_quantiles_df, lineage_colors, segment_order)
    ax_corea.text(-0.02, 1.02, 'C', transform=ax_corea.transAxes, 
                 fontsize=16, fontweight='bold', va='bottom', ha='right')
    
    # Panel D: NP tree (bottom right, first column)
    ax_tree_c = fig.add_subplot(gs[1:3, 1])
    ll_np = load_and_process_tree(path, "NP")
    if ll_np:
        draw_segment_tree(ax_tree_c, ll_np, mrsi_dec, lineage_colors, "NP")
    ax_tree_c.text(-0.02, 1.02, 'D', transform=ax_tree_c.transAxes, 
                  fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_tree_c.set_title('NP', fontsize=12, fontweight='bold', pad=5)
    
    # Panel E: PB2 tree (bottom right, second column)
    ax_tree_d = fig.add_subplot(gs[1:3, 2])
    ll_pb2 = load_and_process_tree(path, "PB2")
    if ll_pb2:
        draw_segment_tree(ax_tree_d, ll_pb2, mrsi_dec, lineage_colors, "PB2")
    ax_tree_d.text(-0.02, 1.02, 'E', transform=ax_tree_d.transAxes, 
                   fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_tree_d.set_title('PB2', fontsize=12, fontweight='bold', pad=5)
    
    # constrained_layout=True in figure creation automatically adjusts subplot positions
    # to ensure all labels, titles, and ticks fit within their allocated space
    return fig


def create_all_segments_figure(path, mrsi_dec, segment_order, lineage_colors):
    """Create figure with all segment trees"""
    # Exclude HA (first segment)
    segments_to_plot = segment_order[1:]
    
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()
    
    for idx, segment in enumerate(segments_to_plot):
        ax = axes[idx]
        ll = load_and_process_tree(path, segment)
        if ll:
            draw_segment_tree(ax, ll, mrsi_dec, lineage_colors, segment)
        ax.set_title(segment, fontweight='bold')
    
    # Add legend in the last subplot
    ax_legend = axes[-1]
    ax_legend.axis('off')
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=lineage_colors['HPAI'], label='HPAI'),
                      Patch(facecolor=lineage_colors['LPAI'], label='LPAI'),
                      Patch(facecolor=lineage_colors['unknown'], label='unknown')]
    ax_legend.legend(handles=legend_elements, loc='center', fontsize=12)
    
    plt.tight_layout()
    return fig


def main(force=False, ind=False):
    """Main function to run the complete analysis"""
    setup_matplotlib()
    
    (clades, rate_shifts, mrsi, mrsi_hpai, mrsi_lpai, segment_order,
     methods_colors, species_colors, clade_colors, lineage_colors, path) = define_constants()
    
    # Convert mrsi to decimal year
    mrsi_dec = (mrsi - datetime(1970, 1, 1)).days / 365.25 + 1970
    
    # Run BEAST commands if needed
    run_beast_commands(path, force=force)
    
    # Load log file
    log_file_path = os.path.join(path, 'combined/HLHxNx.skygrowth.log')
    if not os.path.exists(log_file_path):
        print(f"Error: Log file not found: {log_file_path}")
        return
    
    log_file = pd.read_csv(log_file_path, sep='\t')
    
    # Calculate reassortment rates
    valid_times, all_quantiles, quantiles_to_plot = calculate_reassortment_rates(log_file, mrsi, path)
    
    # Load clade events
    clades_tsv_path = os.path.join(path, 'combined/HLHxNx.skygrowth.clades.tsv')
    if not os.path.exists(clades_tsv_path):
        print(f"Error: Clades TSV file not found: {clades_tsv_path}")
        return
    
    hpai_events = pd.read_csv(clades_tsv_path, sep='\t')
    hpai_events['Time'] = mrsi - pd.to_timedelta(hpai_events['Height'] * 365.25, unit='D')
    hpai_events = hpai_events[hpai_events['Time'] >= datetime(2020, 1, 1)]
    
    # Calculate event distributions
    distr_df, co_rea_df = calculate_event_distribution(hpai_events, segment_order)
    co_rea_quantiles_df = calculate_co_rea_quantiles(co_rea_df)
    
    # Create main figure
    fig_main = create_main_figure(path, mrsi, mrsi_dec, segment_order, lineage_colors, 
                                 species_colors, hpai_events, distr_df, co_rea_df, 
                                 co_rea_quantiles_df, valid_times, all_quantiles, quantiles_to_plot)
    
    # Save main figure
    output_dir = '/Users/nmueller/Documents/github/CoInfection-Material/Figures/'
    os.makedirs(output_dir, exist_ok=True)
    
    if ind:
        fig_main.savefig(os.path.join(output_dir, 'Figure2.pdf'), bbox_inches='tight')
    else:
        fig_main.savefig(os.path.join(output_dir, 'h5n1_reassortment_dependent.pdf'), 
                         bbox_inches='tight')
    
    plt.show()
    plt.close(fig_main)
    
    # Create all segments figure
    fig_segments = create_all_segments_figure(path, mrsi_dec, segment_order, lineage_colors)
    
    if ind:
        fig_segments.savefig(os.path.join(output_dir, 'h5n1_all_segment_trees_skygrowth.pdf'),
                            bbox_inches='tight')
    else:
        fig_segments.savefig(os.path.join(output_dir, 'h5n1_all_segment_trees_dependent.pdf'),
                            bbox_inches='tight')
    
    plt.show()
    plt.close(fig_segments)
    
    print("Analysis complete!")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate H5N1 clade comparison plots')
    parser.add_argument('--force', action='store_true',
                       help='Force rerun of all BEAST commands, even if output files exist')
    parser.add_argument('--ind', action='store_true',
                       help='Use independent model naming for output files')
    
    args = parser.parse_args()
    main(force=args.force, ind=args.ind)

