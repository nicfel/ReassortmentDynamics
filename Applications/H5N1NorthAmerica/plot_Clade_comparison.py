import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
    mpl.rcParams['axes.labelsize'] = 13  # Slightly larger for readability
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
        "HPAI": "#ad5252",            # muted red (Daly City-inspired, matching plot_glm.py)
        "LPAI": "#407499",            # muted blue (Daly City-inspired, matching plot_glm.py)
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
    logcombiner_cmd = f"{beast_path}logcombiner -burnin 20 -log ./out2/HLHxNx.skygrowth.rep*.trees -o ./combined/HLHxNx.skygrowth.trees"
    subprocess.run(logcombiner_cmd, shell=True, cwd=path)
    
    # 2. Summarize network
    print("Running BEAST ReassortmentNetworkSummarize...")
    summarize_cmd = f"{beast_path}applauncher ReassortmentNetworkSummarize -burnin 0 -followSegment 0 -positions MCC ./combined/HLHxNx.skygrowth.trees ./combined/HLHxNx.skygrowth.tree"
    subprocess.run(summarize_cmd, shell=True, cwd=path)
    
    # 3. Combine log files
    print("Running BEAST logcombiner for logs...")
    logcombiner_log_cmd = f"{beast_path}logcombiner -burnin 20 -log ./out2/HLHxNx.skygrowth.rep*.log -o ./combined/HLHxNx.skygrowth.log"
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
    
    # 6. Run ClusterSizeComparison
    cluster_output = os.path.join(path, 'combined/HLHxNx.skygrowth.cluster_comparison.txt')
    if force or not os.path.exists(cluster_output):
        print("Running ClusterSizeComparison...")
        cluster_cmd = f"{beast_path}applauncher ClusterSizeComparison -burnin 0 -clade ./tables/HPAI_LPAI.csv ./combined/HLHxNx.skygrowth.trees ./combined/HLHxNx.skygrowth.cluster_comparison.txt"
        subprocess.run(cluster_cmd, shell=True, cwd=path)
    else:
        print("Cluster size comparison file already exists, skipping...")


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
        
        # Color based on posterior support (gradient matching R code scale but with LPAI blue)
        # R uses: scale_color_gradient2(low="black", mid="#377EB8", high="#377EB8", midpoint=0.7)
        # Using LPAI blue (#407499) instead of green, same scale
        import matplotlib.colors as mcolors
        
        # LPAI blue #407499 is RGB (64, 116, 153)
        lpai_blue_rgb = (64, 116, 153)
        
        # Create gradient: black (0) -> LPAI blue #407499 (0.7) -> LPAI blue #407499 (1)
        if loc <= 0.7:
            # Interpolate from black (0,0,0) to LPAI blue (64, 116, 153) for values 0-0.7
            t = loc / 0.7  # Normalize to 0-1
            r = int(t * lpai_blue_rgb[0])
            g = int(t * lpai_blue_rgb[1])
            b = int(t * lpai_blue_rgb[2])
            col = f'#{r:02x}{g:02x}{b:02x}'
        else:
            # Stay LPAI blue for values > 0.7
            col = '#407499'
        
        # Make linewidth wider for branches with support (scale from 0.5 to 1.0 based on support)
        # Higher support = wider line
        branch_lw = 0.5 + loc*2  # Range from 0.5 to 1.0
        
        # Draw vertical branch (time direction)
        if not isinstance(k, bt.reticulation) and k.parent:
            ax.plot([x, xp], [y, y], color=col, linewidth=branch_lw, zorder=2,
                   solid_capstyle='round', solid_joinstyle='round')
        
        # Draw horizontal node connections
        if k.branchType == 'node' and len(k.children) >= 2:
            left = k.children[-1].y
            right = k.children[0].y
            ax.plot([x, x], [left, right], color=col, linewidth=branch_lw, zorder=2,
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
    # Both legends use x=0.02 for proper horizontal alignment
    legend_x_offset = 0.02
    legend_y_offset = 0.02
    
    legend1 = ax.legend(handles=legend_elements, loc='lower left', fontsize=10, 
                        frameon=False, title='cow isolate', title_fontsize=10,
                        bbox_to_anchor=(legend_x_offset, legend_y_offset))
    
    # Add colorbar for reassortment support (gradient matching R code scale but with LPAI blue)
    # R uses: scale_color_gradient2(low="black", mid="#377EB8", high="#377EB8", midpoint=0.7)
    # Using LPAI blue (#407499) instead of green, same scale
    colors_list = ['black', '#407499']
    n_bins = 100
    cmap = LinearSegmentedColormap.from_list('reassortment', colors_list, N=n_bins)
    
    # Create a scalar mappable for the colorbar
    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    
    # Add colorbar inside the plot area - create axes manually for consistent PDF output
    # Position it higher (y=0.32) to avoid overlap with cow legend at y=0.02
    # Align horizontally with cow legend using the same x offset
    # Get the axes position in figure coordinates
    fig = ax.figure
    ax_pos = ax.get_position()
    
    # Calculate colorbar position in figure coordinates
    # Position in axes coordinates: x=0.01 (slightly left of legend), y=0.32, width=0.15, height=0.03
    cbar_width = 0.15 * ax_pos.width
    cbar_height = 0.03 * ax_pos.height
    cbar_x = ax_pos.x0 + (legend_x_offset - 0.17) * ax_pos.width  # Slightly left of legend
    cbar_y = ax_pos.y0 + 0.32 * ax_pos.height
    
    # Create colorbar axes in figure coordinates for consistent PDF rendering
    cax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
    cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label('posterior support for\nHPAI LPAI reassortment', 
                   fontsize=10, labelpad=8)
    cbar.set_ticks([0, 0.5, 1.0])
    cbar.set_ticklabels(['0', '0.5', '1+'])  # Matching R code breaks and labels
    cbar.ax.tick_params(labelsize=10)
    # Remove background from colorbar axes to ensure proper integration
    cbar.ax.patch.set_facecolor('none')
    cbar.ax.patch.set_edgecolor('none')



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


def plot_reassortment_rates_skygrowth(ax, valid_times, all_quantiles, quantiles_to_plot, mrsi, fromval, toval, lineage_colors):
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
                       alpha=alpha, color=lineage_colors['LPAI'],
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
                   alpha=median_alpha, color=lineage_colors['LPAI'],
                   label='Median', zorder=3, linewidth=0)
    
    ax.set_xlim(fromval, toval)
    ax.set_xlabel('')
    ax.set_ylabel('reassortment rate', labelpad=4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='both')  # Subtle grid
    
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


def calculate_total_quantiles(hpai_events, segment_order):
    """Calculate quantiles for total event counts (non-empty events only) across samples"""
    # Filter for HPAI lineage events only
    hpai_lineage_events = hpai_events[hpai_events['Lineage'] == 'HPAI'].copy()
    
    if len(hpai_lineage_events) == 0:
        return None, None
    
    # Helper function to check if segments are non-empty
    def is_non_empty(seg_str):
        if pd.isna(seg_str):
            return False
        seg_str = str(seg_str).strip('{} ')
        return bool(seg_str)
    
    # Calculate totals per sample
    total_counts = {}
    
    for sample in hpai_lineage_events['Sample'].unique():
        sample_data = hpai_lineage_events[hpai_lineage_events['Sample'] == sample]
        
        # Count unique events within HPAI (Event == "HPAI") with non-empty segments
        hpai_within = sample_data[sample_data['Event'] == 'HPAI']
        hpai_within_non_empty = hpai_within[hpai_within['Segments'].apply(is_non_empty)]
        total_hpai = len(hpai_within_non_empty)
        
        # Count unique events into HPAI from LPAI (Event == "HPAI+LPAI" or "LPAI") with non-empty segments
        hpai_from_lpai = sample_data[sample_data['Event'].isin(['HPAI+LPAI', 'LPAI'])]
        hpai_from_lpai_non_empty = hpai_from_lpai[hpai_from_lpai['Segments'].apply(is_non_empty)]
        total_lpai = len(hpai_from_lpai_non_empty)
        
        if sample not in total_counts:
            total_counts[sample] = {}
        total_counts[sample]['HPAI'] = total_hpai
        total_counts[sample]['LPAI'] = total_lpai
    
    # Calculate quantiles for each event type
    hpai_totals = [total_counts[s]['HPAI'] for s in total_counts.keys()]
    lpai_totals = [total_counts[s]['LPAI'] for s in total_counts.keys()]
    
    if len(hpai_totals) == 0:
        return None, None
    
    q = 0.05
    hpai_lower = np.quantile(hpai_totals, q/2)
    hpai_upper = np.quantile(hpai_totals, 1 - q/2)
    hpai_mean = np.mean(hpai_totals)
    
    lpai_lower = np.quantile(lpai_totals, q/2)
    lpai_upper = np.quantile(lpai_totals, 1 - q/2)
    lpai_mean = np.mean(lpai_totals)
    
    hpai_quantiles = {'mean': hpai_mean, 'lower': hpai_lower, 'upper': hpai_upper}
    lpai_quantiles = {'mean': lpai_mean, 'lower': lpai_lower, 'upper': lpai_upper}
    
    return hpai_quantiles, lpai_quantiles


def create_cluster_size_comparison_plot(ax, path, lineage_colors):
    """Create plot showing probability that children with reassortment > children without, across posterior
    Filtered to resultingClade == HPAI, with separate distributions for incomingClade == LPAI and HPAI"""
    cluster_file_path = os.path.join(path, 'combined/HLHxNx.skygrowth.cluster_comparison.txt')
    
    try:
        cluster_data = pd.read_csv(cluster_file_path, sep='\t')
        
        # Filter for resultingClade == HPAI
        hpai_data = cluster_data[cluster_data['resultingClade'] == 'HPAI'].copy()
        
        if len(hpai_data) == 0:
            ax.text(0.5, 0.5, 'No HPAI resulting clade events found', 
                    transform=ax.transAxes, ha='center', va='center')
            return
        
        # Separate by incomingClade
        lpai_incoming = hpai_data[hpai_data['incomingClade'] == 'LPAI'].copy()
        hpai_incoming = hpai_data[hpai_data['incomingClade'] == 'HPAI'].copy()
        
        # Calculate probabilities for LPAI incoming with bootstrap sampling
        lpai_probs = []
        for iteration in lpai_incoming['iteration'].unique():
            iter_data = lpai_incoming[lpai_incoming['iteration'] == iteration]
            if len(iter_data) > 0:
                # Bootstrap sample (with replacement) from events in this iteration
                n_samples = len(iter_data)
                bootstrap_indices = np.random.choice(iter_data.index, size=n_samples, replace=True)
                bootstrap_data = iter_data.loc[bootstrap_indices]
                
                # Exclude ties (where leafsWith == leafsWithout) from calculation
                non_tie_data = bootstrap_data[bootstrap_data['leafsWith'] != bootstrap_data['leafsWithout']]
                
                if len(non_tie_data) > 0:
                    # Count strict greater cases in bootstrapped sample (excluding ties)
                    with_greater = (non_tie_data['leafsWith'] > non_tie_data['leafsWithout']).sum()
                    total_comparisons = len(non_tie_data)
                    prob_with_greater = with_greater / total_comparisons if total_comparisons > 0 else 0
                    lpai_probs.append(prob_with_greater)
        
        # Calculate probabilities for HPAI incoming with bootstrap sampling
        hpai_probs = []
        for iteration in hpai_incoming['iteration'].unique():
            iter_data = hpai_incoming[hpai_incoming['iteration'] == iteration]
            if len(iter_data) > 0:
                # Bootstrap sample (with replacement) from events in this iteration
                n_samples = len(iter_data)
                bootstrap_indices = np.random.choice(iter_data.index, size=n_samples, replace=True)
                bootstrap_data = iter_data.loc[bootstrap_indices]
                
                # Exclude ties (where leafsWith == leafsWithout) from calculation
                non_tie_data = bootstrap_data[bootstrap_data['leafsWith'] != bootstrap_data['leafsWithout']]
                
                if len(non_tie_data) > 0:
                    # Count strict greater cases in bootstrapped sample (excluding ties)
                    with_greater = (non_tie_data['leafsWith'] > non_tie_data['leafsWithout']).sum()
                    total_comparisons = len(non_tie_data)
                    prob_with_greater = with_greater / total_comparisons if total_comparisons > 0 else 0
                    hpai_probs.append(prob_with_greater)
        
        # Prepare data for violin plot
        plot_data = []
        plot_labels = []
        
        if len(lpai_probs) > 0:
            plot_data.append(lpai_probs)
            plot_labels.append('LPAI')
        
        if len(hpai_probs) > 0:
            plot_data.append(hpai_probs)
            plot_labels.append('HPAI')
        
        # Create violin plot
        if len(plot_data) > 0:
            positions = range(len(plot_data))
            parts = ax.violinplot(plot_data, positions=positions, widths=0.7, 
                                 showmeans=False, showmedians=False, showextrema=False)
            
            # Color the violins
            for i, pc in enumerate(parts['bodies']):
                pc.set_facecolor(lineage_colors[plot_labels[i]])
                pc.set_alpha(0.7)
                pc.set_edgecolor('black')
                pc.set_linewidth(1)
            
            # Style the quartile lines
            if 'cquantiles' in parts:
                parts['cquantiles'].set_color('black')
                parts['cquantiles'].set_linewidth(1)
            
            # Set x-axis labels
            ax.set_xticks(positions)
            ax.set_xticklabels(plot_labels)
        
        # Add dashed horizontal line at y=0.5 - more visible reference line
        ax.axhline(y=0.5, color='#666666', linestyle='--', linewidth=2, alpha=0.8, zorder=0)

        ax.set_xlabel('incoming', labelpad=3)
        ax.set_ylabel('P(reassortant > non-reassortant)', labelpad=4)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='y')  # More subtle grid
        ax.set_ylim(0, 1)
        
    except FileNotFoundError:
        ax.text(0.5, 0.5, 'Cluster comparison file not found', 
                transform=ax.transAxes, ha='center', va='center')
    except Exception as e:
        ax.text(0.5, 0.5, f'Error loading cluster data:\n{str(e)}', 
                transform=ax.transAxes, ha='center', va='center', fontsize=8)
        import traceback
        traceback.print_exc()


def plot_co_reassortment(ax, co_rea_quantiles_df, lineage_colors, segment_order, hpai_events=None):
    """Plot co-reassortment events by segment"""
    # Filter for HPAI lineage only
    data = co_rea_quantiles_df[co_rea_quantiles_df['Lineage'] == 'HPAI'].copy()
    
    if len(data) == 0:
        return
    
    desired_order = segment_order[1:]  # skip HA
    segments = [seg for seg in desired_order if seg in set(data['segment'])]
    
    # Add "Total" to segments if hpai_events is provided
    if hpai_events is not None:
        segments = segments + ['Total']
    
    x_pos = np.arange(len(segments))
    width = 0.35
    dodge_offset = -0.3  # Match R's position_dodge(-0.3)
    
    # Pre-calculate total quantiles if hpai_events is provided
    total_quantiles = None
    if hpai_events is not None:
        total_quantiles = calculate_total_quantiles(hpai_events, segment_order)
    
    for evname in ['HPAI', 'LPAI']:
        ev_data = data[data['evname'] == evname]
        if len(ev_data) == 0:
            continue
        
        means = []
        lowers = []
        uppers = []
        
        # Process individual segments
        for s in segments:
            if s == 'Total':
                # Use pre-calculated totals
                if total_quantiles is not None:
                    hpai_quantiles, lpai_quantiles = total_quantiles
                    if hpai_quantiles is not None and lpai_quantiles is not None:
                        if evname == 'HPAI':
                            means.append(hpai_quantiles['mean'])
                            lowers.append(hpai_quantiles['lower'])
                            uppers.append(hpai_quantiles['upper'])
                        else:  # LPAI
                            means.append(lpai_quantiles['mean'])
                            lowers.append(lpai_quantiles['lower'])
                            uppers.append(lpai_quantiles['upper'])
                    else:
                        means.append(0)
                        lowers.append(0)
                        uppers.append(0)
                else:
                    means.append(0)
                    lowers.append(0)
                    uppers.append(0)
            else:
                # Regular segment processing
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
    ax.set_ylabel('events', labelpad=3)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(segments)
    ax.legend(loc='upper left', fontsize=10, title='segment origin', frameon=False,
              bbox_to_anchor=(0, 1.15))  # Slightly lower to avoid overlap
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.2, linestyle='-', linewidth=0.5, axis='y')  # More subtle grid


def create_main_figure(path, mrsi, mrsi_dec, segment_order, lineage_colors, species_colors, 
                       hpai_events, distr_df, co_rea_df, co_rea_quantiles_df, 
                       valid_times, all_quantiles, quantiles_to_plot):
    """Create the main composite figure (Figure 4)"""
    fig = plt.figure(figsize=(14, 10.5), constrained_layout=True)
    # Add padding to ensure panel labels outside plot area are visible
    # Increased padding to accommodate labels positioned at (-0.02, 1.02)
    fig.set_constrained_layout_pads(w_pad=0.08, h_pad=0.08)
    # Layout: Compact and balanced - 3x3 grid
    # A: HA tree (left column, spans rows 0-1) - main focal point, largest panel
    # B: Reassortment rates (top right, row 0, spans cols 1-2) - time series needs width
    # C: Co-reassortment (bottom left, row 2, col 0) - bar plot needs width for categories
    # D: NP tree (middle right, row 1, col 1) - compact tree
    # E: PB2 tree (middle right, row 1, col 2) - compact tree  
    # F: Cluster comparison (bottom right, row 2, spans cols 1-2) - violin plot, compact width
    # Use constrained_layout to ensure all content fits within subplot boundaries
    # Width ratios: C (col 0) gets more space than F (cols 1-2 combined)
    gs = gridspec.GridSpec(3, 3, width_ratios=[4, 1, 1], height_ratios=[1.8, 3.2, 1.4], 
                          figure=fig)

    # Consistent panel label style - positioned OUTSIDE plot area (standard scientific figure style)
    # Labels positioned to the left and above each subplot using transform=transAxes
    # Position: (-0.02, 1.02) means 2% to left of left edge, 2% above top edge
    # Alignment: va='bottom', ha='right' aligns label's bottom-right corner at that position
    panel_label_kwargs = {'fontsize': 16, 'fontweight': 'bold', 'va': 'bottom', 'ha': 'right'}
    panel_label_pos = (-0.02, 1.02)  # Standard: outside top-left corner
    
    # Panel A: HA tree with clades (left column, spans rows 0-1 - main focal point)
    ax_tree = fig.add_subplot(gs[0:2, 0])
    ll_ha = load_and_process_tree(path, "HA")
    if ll_ha:
        draw_ha_tree_with_clades(ax_tree, ll_ha, mrsi, mrsi_dec, lineage_colors, species_colors)
    ax_tree.text(panel_label_pos[0], panel_label_pos[1], 'A', transform=ax_tree.transAxes, **panel_label_kwargs)
    
    # Panel B: Reassortment rates (top right, row 0, spans cols 1-2)
    ax_rates = fig.add_subplot(gs[0, 1:])
    if valid_times is not None and len(valid_times) > 0:
        fromval = mrsi_dec - 10  # Adjust as needed
        toval = mrsi_dec
        plot_reassortment_rates_skygrowth(ax_rates, valid_times, all_quantiles, quantiles_to_plot, 
                                         mrsi, fromval, toval, lineage_colors)
    ax_rates.text(panel_label_pos[0], panel_label_pos[1], 'B', transform=ax_rates.transAxes, **panel_label_kwargs)
    
    # Panel C: Co-reassortment plot (bottom left, row 2, col 0)
    ax_corea = fig.add_subplot(gs[2, 0])
    plot_co_reassortment(ax_corea, co_rea_quantiles_df, lineage_colors, segment_order, hpai_events)
    ax_corea.text(panel_label_pos[0], panel_label_pos[1], 'C', transform=ax_corea.transAxes, **panel_label_kwargs)
    
    # Panel D: NP tree (middle right, row 1, col 1)
    ax_tree_c = fig.add_subplot(gs[1, 1])
    ll_np = load_and_process_tree(path, "NP")
    if ll_np:
        draw_segment_tree(ax_tree_c, ll_np, mrsi_dec, lineage_colors, "NP")
    ax_tree_c.text(panel_label_pos[0], panel_label_pos[1], 'D', transform=ax_tree_c.transAxes, **panel_label_kwargs)
    
    # Panel E: PB2 tree (middle right, row 1, col 2)
    ax_tree_d = fig.add_subplot(gs[1, 2])
    ll_pb2 = load_and_process_tree(path, "PB2")
    if ll_pb2:
        draw_segment_tree(ax_tree_d, ll_pb2, mrsi_dec, lineage_colors, "PB2")
    ax_tree_d.text(panel_label_pos[0], panel_label_pos[1], 'E', transform=ax_tree_d.transAxes, **panel_label_kwargs)
    
    # Panel F: Cluster size comparison (bottom right, row 2, spans cols 1-2 - compact histogram)
    ax_cluster = fig.add_subplot(gs[2, 1:])
    create_cluster_size_comparison_plot(ax_cluster, path, lineage_colors)
    ax_cluster.text(panel_label_pos[0], panel_label_pos[1], 'F', transform=ax_cluster.transAxes, **panel_label_kwargs)
    
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
        fig_main.savefig(os.path.join(output_dir, 'Figure4.pdf'), bbox_inches='tight')
    else:
        fig_main.savefig(os.path.join(output_dir, 'Figure4.pdf'), 
                         bbox_inches='tight')
    
    # Don't show figures - just save them
    plt.close(fig_main)
    
    # Create all segments figure
    fig_segments = create_all_segments_figure(path, mrsi_dec, segment_order, lineage_colors)
    
    fig_segments.savefig(os.path.join(output_dir, 'h5n1_all_segment_trees_skygrowth.pdf'),
                            bbox_inches='tight')
    
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

