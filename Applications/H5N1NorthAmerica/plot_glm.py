import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
from scipy.stats import gaussian_kde
import os
import random
from itertools import permutations
import sys
import pandas as pd
import subprocess
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
sys.path.append('/Users/nmueller/Documents/github/CoInfection-Material/Applications/NetworkViz/')
import baltic_bacter as bt


def setup_matplotlib():
    """Set up matplotlib defaults"""
    typeface = 'helvetica'
    mpl.rcParams['font.weight'] = 300
    mpl.rcParams['axes.labelweight'] = 300
    mpl.rcParams['font.family'] = typeface
    mpl.rcParams['font.size'] = 12  # Consistent base font size
    mpl.rcParams['axes.labelsize'] = 12
    mpl.rcParams['axes.titlesize'] = 12
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['legend.fontsize'] = 10


def define_constants():
    """Define all constants used in the analysis"""
    colours = {'red', 'green'}
    segments = ['HA', 'NA', 'MP', 'NS', 'NP', 'PB1', 'PB2', 'PA']
    new_order = [0, 1, 2, 4, 7, 5, 6, 3]
    # Daly City-inspired color scheme: darker, muted tones
    colour_cycle = ['#5B7A9F', '#8B7AA6', '#6B9568', '#C89878', '#9B8874', '#6B8583']
    path = '/Users/nmueller/Documents/github/CoInfection-Material/Applications/H5N1NorthAmerica/'
    linewidth = 0.5
    approach = 'glm'

    return colours, segments, new_order, colour_cycle, path, linewidth, approach


def run_beast_commands(path, force=False):
    """Run BEAST commands to process the data"""
    # Define BEAST path
    beast_path = "/Applications/BEAST\\ 2.7.7/bin/"
    
    # Create combined directory if it doesn't exist
    combined_dir = os.path.join(path, 'combined')
    os.makedirs(combined_dir, exist_ok=True)
    
    # Check if main output files exist and skip if not forcing
    trees_output = os.path.join(path, 'combined/HPAI_HLHxNx.glm.trees')
    tree_output = os.path.join(path, 'combined/HPAI_HLHxNx.glm.tree')
    log_output = os.path.join(path, 'combined/HPAI_HLHxNx.glm.log')
    b113_output = os.path.join(path, 'combined/b113_glm.tsv')
    d11_output = os.path.join(path, 'combined/d11_glm.tsv')
    
    if not force and all(os.path.exists(f) for f in [trees_output, tree_output, log_output, b113_output, d11_output]):
        print("Main BEAST output files already exist. Skipping main BEAST commands.")
        print("Use force=True to rerun all commands.")
        # Don't return here - continue to check segment and cluster processing
    else:
        print("Some main BEAST files missing, will process as needed...")
    
    if force:
        print("Force mode enabled - rerunning all BEAST commands...")
    
    # 1. Combine tree files
    if force or not os.path.exists(trees_output):
        print("Running BEAST logcombiner for trees...")
        logcombiner_cmd = f"{beast_path}logcombiner -burnin 20 -log ./out/HPAI_HLHxNx.glm.rep*.trees -o ./combined/HPAI_HLHxNx.glm.trees"
        subprocess.run(logcombiner_cmd, shell=True, cwd=path)
    else:
        print("Tree combination file already exists, skipping...")
    
    # 2. Summarize network
    if force or not os.path.exists(tree_output):
        print("Running BEAST ReassortmentNetworkSummarize...")
        summarize_cmd = f"{beast_path}applauncher ReassortmentNetworkSummarize -burnin 0 -followSegment 0 -positions MCC ./combined/HPAI_HLHxNx.glm.trees ./combined/HPAI_HLHxNx.glm.tree"
        subprocess.run(summarize_cmd, shell=True, cwd=path)
    else:
        print("Network summary file already exists, skipping...")
    
    # 3. Combine log files
    if force or not os.path.exists(log_output):
        print("Running BEAST logcombiner for logs...")
        logcombiner_log_cmd = f"{beast_path}logcombiner -burnin 20 -log ./out/HPAI_HLHxNx.glm.rep*.log -o ./combined/HPAI_HLHxNx.glm.log"
        subprocess.run(logcombiner_log_cmd, shell=True, cwd=path)
    else:
        print("Log combination file already exists, skipping...")
    
    # 4. Extract clade heights
    if force or not os.path.exists(b113_output):
        print("Extracting clade heights for B3.13...")
        clade_heights_cmd1 = f"{beast_path}applauncher GetCladeHeightsFromNetwork -burnin 0 -tree ./combined/HPAI_HLHxNx.glm.trees -clade ./tables/cow_clade.csv -out ./combined/b113_glm.tsv"
        subprocess.run(clade_heights_cmd1, shell=True, cwd=path)
    else:
        print("B3.13 clade heights file already exists, skipping...")
    
    if force or not os.path.exists(d11_output):
        print("Extracting clade heights for D1.1...")
        clade_heights_cmd2 = f"{beast_path}applauncher GetCladeHeightsFromNetwork -burnin 0 -tree ./combined/HPAI_HLHxNx.glm.trees -clade ./tables/d11.csv -out ./combined/d11_glm.tsv"
        subprocess.run(clade_heights_cmd2, shell=True, cwd=path)
    else:
        print("D1.1 clade heights file already exists, skipping...")
    
    # 5. Process segments
    print("Processing segments...")
    segment_order = ["HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA"]
    for s, segment in enumerate(segment_order):
        segment_output = os.path.join(path, f'combined/HPAI_HLHxNx.glm.{segment}.trees')
        if force or not os.path.exists(segment_output):
            print(f"Processing segment {segment}...")
            segment_cmd = f"{beast_path}applauncher MarkCladesFromCladeFile -burnin 0 -followSegment {s} -printSegment {s} -tree ./combined/HPAI_HLHxNx.glm.trees -clade ./tables/HPAI_LPAI.csv -out ./combined/HPAI_HLHxNx.glm.{segment}.trees"
            subprocess.run(segment_cmd, shell=True, cwd=path)
        else:
            print(f"Segment {segment} file already exists, skipping...")
    
    # 6. Run ClusterSizeComparison
    cluster_output = os.path.join(path, 'combined/HPAI_HLHxNx.glm.cluster_comparison.txt')
    if force or not os.path.exists(cluster_output):
        print("Running ClusterSizeComparison...")
        cluster_cmd = f"{beast_path}applauncher ClusterSizeComparison -burnin 0 ./combined/HPAI_HLHxNx.glm.trees ./combined/HPAI_HLHxNx.glm.cluster_comparison.txt"
        subprocess.run(cluster_cmd, shell=True, cwd=path)
    else:
        print("Cluster size comparison file already exists, skipping...")


def calculate_clade_probabilities(path, force=False):
    """Calculate reassortment event probabilities for clades"""
    output_path = os.path.join(path, 'combined/HPAI_HLHxNx.glm.cladeprobs.csv')
    
    # Check if output file already exists
    if not force and os.path.exists(output_path):
        print("Clade probabilities file already exists. Loading existing data...")
        return pd.read_csv(output_path)
    
    if force:
        print("Force mode enabled - recalculating clade probabilities...")
    
    # Define constants
    clades = ["B3.13", "D1.1"]
    rate_shift_str = '0 0.105936073059361 0.211872146118721 0.317808219178082 0.423744292237443 0.529680365296804 0.635616438356164 0.741552511415525 0.847488584474886 0.953424657534247 1.05936073059361 1.16529680365297 1.27123287671233 1.37716894977169 1.48310502283105 1.58904109589041 1.69497716894977 1.80091324200913 1.90684931506849 2.01278538812785 2.11872146118721 2.22465753424658 2.33059360730594 2.4365296803653 2.54246575342466 2.64840182648402 2.75433789954338 2.86027397260274 2.9662100456621 3.07214611872146 3.17808219178082 3.28401826484018 3.38995433789954 3.4958904109589 3.60182648401827 3.70776255707763 3.81369863013699 3.91963470319635 4.02557077625571 4.13150684931507'
    rate_shifts = np.array([float(x) for x in rate_shift_str.split()])
    
    print("Reading clade heights and log data...")
    # Read clade heights
    clade_cow_heights = pd.read_csv(os.path.join(path, 'combined/b113_glm.tsv'), sep='\t')
    clade_d11_heights = pd.read_csv(os.path.join(path, 'combined/d11_glm.tsv'), sep='\t')
    log_file = pd.read_csv(os.path.join(path, 'combined/HPAI_HLHxNx.glm.log'), sep='\t')
    
    data = []
    
    print("Calculating probabilities for each clade...")
    for cl in clades:
        if cl == "B3.13":
            clade_heights = clade_cow_heights
        elif cl == "D1.1":
            clade_heights = clade_d11_heights
        
        # Loop over the posterior
        for l in range(min(len(clade_heights), len(log_file))):
            # Get the timings of the HA segment
            min_time = clade_heights.iloc[l, 1]  # Second column (index 1)
            max_time = min_time + 0.5
            
            first_interval = len([x for x in rate_shifts if x <= min_time]) - 1
            last_interval = len([x for x in rate_shifts if x <= max_time]) - 1
            
            curr_time = min_time
            weighted = 0.0
            
            for i in range(first_interval, last_interval + 1):
                if i >= len(rate_shifts) - 1:
                    break
                    
                next_time = min(rate_shifts[i + 1], max_time)
                
                # Get the reassortment rates of this interval at the beginning and end
                r_start = log_file.iloc[l, log_file.columns.get_loc(f'reassortmentRate.{i}')]
                r_end = log_file.iloc[l, log_file.columns.get_loc(f'reassortmentRate.{i + 1}')]
                
                # Calculate the growth rate for this interval
                growth = (r_start - r_end) / (rate_shifts[i + 1] - rate_shifts[i])
                
                timediff1 = curr_time - rate_shifts[i]
                timediff2 = next_time - rate_shifts[i]
                
                if growth == 0.0:
                    weighted = weighted + (next_time - curr_time) * np.exp(r_start)
                else:
                    weighted = weighted + np.exp(r_start) / (-growth) * (
                        np.exp(-growth * timediff2) - np.exp(-growth * timediff1)
                    )
                curr_time = next_time
            
            # Compute probability of no event over interval
            no_event_prob = 1 - (1 - 0.5**(8-1)) * np.exp(-weighted)
            
            data.append({
                'no_event_prob': no_event_prob,
                'min_time': min_time,
                'max_time': max_time,
                'mean_rate': weighted / (max_time - min_time),
                'weighted': weighted,
                'clade': cl
            })
    
    # Convert to DataFrame and save
    data_df = pd.DataFrame(data)
    data_df.to_csv(output_path, index=False)
    print(f"Saved clade probabilities to {output_path}")
    
    return data_df


def process_glm_data(path, force=False):
    """Main function to process GLM data - replaces convert_GLM_data.R"""
    print("Starting GLM data processing...")
    
    # Run BEAST commands
    run_beast_commands(path, force=force)
    
    # Calculate clade probabilities
    clade_data = calculate_clade_probabilities(path, force=force)
    
    print("GLM data processing complete!")
    return clade_data


def load_and_process_tree(bt, path, approach):
    """Load and process the phylogenetic tree"""
    file = 'combined/HPAI_HLHxNx.' + approach + '.tree'
    tree_path = os.path.join(path, file)
    print(tree_path)
    ll = bt.loadNexus(tree_path, date_fmt='%Y-%m-%d', verbose=False)
    ll.drawTree()
    return ll

def setup_tree_plot(ax, ll):
    """Configure the tree plot axes and timeline"""
    ax.set_facecolor('w')
    [ax.spines[loc].set_visible(False) for loc in ax.spines if loc != 'left']

    fromval = float(int(ll.root.absoluteTime + ll.treeHeight)) - 6
    toval = float(int(ll.root.absoluteTime + ll.treeHeight)) + 0.5
    timewidth = 0.5

    for i in np.arange(fromval, toval, 2 * timewidth):
        ax.axhspan(i, i + timewidth, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)

    return fromval, toval, timewidth

def initialize_traits(ll):
    """Initialize traits for all tree objects"""
    for k in ll.Objects:
        k.traits['re'] = 0

def assign_reassortment_colors(ll, colour_cycle):
    """Assign colors to reassortment events"""
    curr_traits_number = 0

    for k in sorted(ll.Objects, key=lambda w: w.height):
        if hasattr(k, 'contribution'):
            random_number = random.randint(0, len(colour_cycle) - 1)
            while random_number == k.traits['re']:
                random_number = random.randint(0, len(colour_cycle) - 1)

            subtree = ll.traverse_tree(k.children[-1], include_condition=lambda w: True)
            for w in subtree:
                w.traits['re'] = random_number

def draw_tree_branches(ax, ll, colour_cycle, linewidth, fromval):
    """Draw all tree branches and nodes with improved aesthetics"""
    reassortment_events = []

    # Define color palette - using consistent colors for all lines
    main_branch_color = '#2c3e50'  # Dark blue-grey for main branches
    secondary_branch_color = '#95a5a6'  # Light grey for secondary branches
    reassortment_color = '#e74c3c'  # Red for reassortment events
    leaf_outer_color = '#34495e'  # Dark grey for leaf outer
    leaf_inner_color = '#ecf0f1'  # Light grey for leaf inner

    for k in ll.Objects:
        x = k.absoluteTime
        xp = k.parent.absoluteTime
        if xp != None:
            xp = max(xp, fromval + 0.000001)
        y = k.y
        col = colour_cycle[k.traits['re'] % len(colour_cycle)]

        if isinstance(k, bt.reticulation) == False:
            col_lin = secondary_branch_color
            lw_scale = 1.2
            alpha = 0.7
            if k.traits['seg0'] == 'true':
                col_lin = col
                lw_scale = 2.5
                alpha = 1.0

            # Swap: horizontal branches become vertical
            ax.plot([y, y], [x, xp], color=col_lin, lw=linewidth * lw_scale,
                   solid_capstyle='round', solid_joinstyle='round', alpha=1, zorder=2)
        else:
            # Swap: horizontal branches become vertical - reassortment branches use main color
            ax.plot([y, y], [x, xp], color=col, lw=linewidth * 1.5,
                   ls='--', solid_capstyle='round', solid_joinstyle='round',
                   alpha=1, zorder=1)

        if k.branchType == 'node':
            left, right = k.children[-1].y, k.children[0].y

            col_lin1 = secondary_branch_color
            lw_scale1 = 1.2
            alpha1 = 0.7
            col_lin2 = secondary_branch_color
            lw_scale2 = 1.2
            alpha2 = 0.7

            # if k.children[-1].traits['seg0'] == 'true':
            col_lin1 = col
            lw_scale1 = 2.5
            alpha1 = 1.0

            # if k.children[0].traits['seg0'] == 'true':
            col_lin2 = col
            lw_scale2 = 2.5
            alpha2 = 1.0

            # Swap: vertical lines become horizontal with rounded joins
            ax.plot([left, k.y], [x, x], color=col_lin1, lw=linewidth * lw_scale1,
                   solid_capstyle='round', solid_joinstyle='round', alpha=alpha1, zorder=2)
            ax.plot([k.y, right], [x, x], color=col_lin2, lw=linewidth * lw_scale2,
                   solid_capstyle='round', solid_joinstyle='round', alpha=alpha2, zorder=2)

        elif isinstance(k, bt.leaf):
            # Swap x and y - improved leaf nodes
            ax.scatter(y, x, s=30, facecolor=leaf_outer_color, edgecolor='none', zorder=4)
            ax.scatter(y, x, s=15, facecolor=col, edgecolor='none', zorder=5)

        elif isinstance(k, bt.reticulation):
            segs = sorted(map(int, k.traits['segments']))

            reassortment_events.append({
                'time': x,
                'segments': segs,
                'posterior': k.traits.get('posterior', 1.0)
            })

            # Swap x and y - smaller reassortment nodes
            ax.scatter(k.target.y, x, s=20, facecolor=reassortment_color,
                      edgecolor='white', linewidth=0.5, zorder=4, alpha=0.9)
            ax.scatter(k.target.y, x, s=6, facecolor=reassortment_color, edgecolor='none', zorder=5)
            ax.plot([y, k.target.y], [x, x], color=col, lw=linewidth * 1.5,
                   ls='-', solid_capstyle='round', solid_joinstyle='round', alpha=1, zorder=1)

            for i in range(len(segs)):
                name = segs[i]
                c = 'black'
                o = 1 / 20.
                posterior_val = round(k.traits['posterior'], 2)

    return reassortment_events

def finalize_tree_plot(ax, ll, fromval, toval):
    """Apply final formatting to the tree plot"""
    ax.set_xticks([])
    ax.set_xlim(ll.ySpan * 1.01, -ll.ySpan * 0.05)  # Inverted x-axis
    ax.set_ylim(toval, fromval)  # toval (recent) at bottom, fromval (past) at top - inverted

    # Add y-axis with time labels
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_linewidth(0.5)
    ax.set_ylabel('Time (years)')

    # Set y-axis ticks
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', labelsize=10)

def load_log_data(path, approach):
    """Load and process the log file data"""
    log_file_path = os.path.join(path, 'combined/HPAI_HLHxNx.' + approach + '.log')
    print(sys.executable)
    print(log_file_path)
    log_data = pd.read_csv(log_file_path, sep='\t')
    return log_data

def calculate_rate_quantiles(log_data, mrsi, path):
    """Calculate quantiles for reassortment rates over time"""
    # Extract rate shifts from XML file
    xml_file_path = os.path.join(path, 'xmls/HPAI_HLHxNx.glm.rep0.xml')
    with open(xml_file_path, 'r') as xml_file:
        xml_content = xml_file.read()
    
    import re
    rate_shifts_match = re.findall(r'<stateNode id="rateShifts" spec="RealParameter" value="([^"]+)"/>', xml_content)
    rate_shifts = np.array([float(x) for x in rate_shifts_match[0].split()])
    
    time_points = [mrsi - shift for shift in rate_shifts]
    
    # Define quantiles for proper confidence intervals
    # 95% CI: 2.5% to 97.5%, 90% CI: 5% to 95%, 80% CI: 10% to 90%, etc.
    quantiles_to_plot = [0.025, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.975]
    
    # Create colors that go from light to dark (light blue to dark blue)
    n_quantiles = len(quantiles_to_plot)
    colors_quantiles = []
    for i in range(n_quantiles):
        # Create gradient from light blue to dark blue
        intensity = 1.0 - (i / (n_quantiles - 1)) * 0.7  # From 1.0 to 0.3
        colors_quantiles.append((0.3, 0.47, 0.65, intensity))  # RGBA format
    
    all_quantiles = []
    all_quantiles_ne = []
    valid_times = []
    
    medians_ne = []
    lower_95_ne = []
    upper_95_ne = []
    lower_50_ne = []
    upper_50_ne = []
    
    for i, shift in enumerate(rate_shifts):
        if i >= len(time_points):
            break
        
        col_name = f'reassortmentRate.{i}'
        col_name2 = f'logNe.{i+1}'
        rates_at_time = log_data[col_name].values
        ne_at_time = log_data[col_name2].values
        
        # Calculate all quantiles for this time point (rates)
        time_quantiles = []
        for q in quantiles_to_plot:
            time_quantiles.append(np.quantile(rates_at_time, q))
        all_quantiles.append(time_quantiles)
        
        # Calculate all quantiles for this time point (Ne)
        time_quantiles_ne = []
        for q in quantiles_to_plot:
            time_quantiles_ne.append(np.quantile(ne_at_time, q))
        all_quantiles_ne.append(time_quantiles_ne)
        
        valid_times.append(time_points[i])
        
        # Keep the original NE calculations for compatibility
        medians_ne.append(np.quantile(ne_at_time, 0.5))
        lower_95_ne.append(np.quantile(ne_at_time, 0.025))
        upper_95_ne.append(np.quantile(ne_at_time, 0.975))
        lower_50_ne.append(np.quantile(ne_at_time, 0.25))
        upper_50_ne.append(np.quantile(ne_at_time, 0.75))
    
    # Convert to numpy arrays for easier handling
    all_quantiles = np.array(all_quantiles)  # Shape: (time_points, quantiles)
    all_quantiles_ne = np.array(all_quantiles_ne)  # Shape: (time_points, quantiles)
    
    # For backward compatibility, also return the original format
    medians = all_quantiles[:, -1]  # 50% quantile (median) - last element
    lower_95 = all_quantiles[:, 0]  # 95% quantile
    upper_95 = all_quantiles[:, 0]  # 95% quantile (same as lower for now)
    lower_50 = all_quantiles[:, -1]  # 50% quantile (median) - last element
    upper_50 = all_quantiles[:, -1]  # 50% quantile (same as lower for now)
    
    return (valid_times, medians, lower_95, upper_95, lower_50, upper_50,
            medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, 
            time_points, all_quantiles, all_quantiles_ne, quantiles_to_plot, colors_quantiles)

def plot_reassortment_rates(ax, valid_times, medians, lower_95, upper_95, lower_50, upper_50, fromval, mrsi, path, all_quantiles=None, quantiles_to_plot=None, colors_quantiles=None, timewidth=0.5):
    """Plot reassortment rates with confidence intervals"""
    offset_y = 500.
    stretch_y = -25
    xml_file_path = os.path.join(path, 'xmls/HPAI_HLHxNx.glm.rep0.xml')

    with open(xml_file_path, 'r') as xml_file:
        xml_content = xml_file.read()

    # Extract actual predictors from the XML (only 3 exist: lpai, hpai, total)
    import re

    lpai_values = re.findall(r'<stateNode id="lpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    hpai_values = re.findall(r'<stateNode id="hpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    total_values = re.findall(r'<stateNode id="total" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)

    lpai_values = np.array([float(x) for x in lpai_values[0].split()])
    hpai_values = np.array([float(x) for x in hpai_values[0].split()])
    total_values = np.array([float(x) for x in total_values[0].split()])

    # set all values 0 to NaN
    lpai_values[lpai_values == 0] = np.nan
    hpai_values[hpai_values == 0] = np.nan
    total_values[total_values == 0] = np.nan

    # standardize the values
    lpai_values = (lpai_values - np.nanmean(lpai_values)) / np.nanstd(lpai_values)
    hpai_values = (hpai_values - np.nanmean(hpai_values)) / np.nanstd(hpai_values)
    total_values = (total_values - np.nanmean(total_values)) / np.nanstd(total_values)

    # Ensure arrays have the same length as valid_times
    target_length = len(valid_times)
    lpai_values = lpai_values[:target_length]
    hpai_values = hpai_values[:target_length]
    total_values = total_values[:target_length]
    
    medians = np.array(medians)
    lower_95 = np.array(lower_95)
    upper_95 = np.array(upper_95)
    lower_50 = np.array(lower_50)
    upper_50 = np.array(upper_50)
    
    # Plot confidence intervals as shaded ribbons
    if all_quantiles is not None and quantiles_to_plot is not None and colors_quantiles is not None:
        # Define confidence intervals and their corresponding quantile pairs
        confidence_intervals = [
            (0.95, 0.025, 0.975),  # 95% CI: 2.5% to 97.5%
            (0.90, 0.05, 0.95),    # 90% CI: 5% to 95%
            (0.80, 0.1, 0.9),      # 80% CI: 10% to 90%
            (0.60, 0.2, 0.8),      # 60% CI: 20% to 80%
            (0.40, 0.3, 0.7),      # 40% CI: 30% to 70%
            (0.20, 0.4, 0.6),      # 20% CI: 40% to 60%
        ]
        
        # Plot each confidence interval as a shaded ribbon
        for ci_level, lower_q, upper_q in confidence_intervals:
            # Find the indices for these quantiles
            lower_idx = quantiles_to_plot.index(lower_q)
            upper_idx = quantiles_to_plot.index(upper_q)
            
            # Get the quantile values
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

        # Plot the median as a thick shaded ribbon
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
    else:
        # Fallback to original plotting if new parameters not provided
        ax.fill_between(valid_times, lower_95, upper_95, alpha=0.2, color='#4E79A7', label='95% CI', zorder=2)
        ax.fill_between(valid_times, lower_50, upper_50, alpha=0.4, color='#4E79A7', label='50% CI', zorder=2)
        ax.plot(valid_times, medians, color='#4E79A7', linewidth=2, label='Median', zorder=2)

    # Plot only LPAI and HPAI predictors with optimal scaling for reassortment rates
    reassortment_scale = 2.97
    reassortment_offset = -2.50
    
    ax.plot(valid_times, lpai_values * reassortment_scale + reassortment_offset, color='#2ca02c', linewidth=1.5, label='LPAI', linestyle='--')
    ax.plot(valid_times, hpai_values * reassortment_scale + reassortment_offset, color='#d62728', linewidth=1.5, label='HPAI', linestyle='--')

    # Add timeline shading to match the tree
    for i in np.arange(fromval, mrsi, 2 * timewidth):
        ax.axvspan(i, i + timewidth, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)

    ax.set_xlim(fromval, mrsi)

    # Set x-axis ticks to only show full years
    year_ticks = np.arange(np.ceil(fromval), np.floor(mrsi) + 1, 1.0)
    ax.set_xticks(year_ticks)
    ax.set_xticklabels([int(year) for year in year_ticks])

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set logarithmic y-axis with custom breaks and labels
    # Define the tick positions (these correspond to log values)
    tick_positions = [np.log(0.025), np.log(0.05), np.log(0.1), np.log(0.2), np.log(0.4), 
                     np.log(0.8), np.log(1.6), np.log(3.2)]
    
    tick_labels = ["0.025","0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "3.2"]
    
    ax.set_yticks(tick_positions)
    ax.set_yticklabels(tick_labels)
    ax.set_ylabel('reassortment rate')
    ax.set_ylim(np.log(0.0125), np.log(6.4))

    return offset_y, stretch_y


def plot_ne(ax_ne, valid_times, medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, fromval, mrsi, path, all_quantiles_ne=None, quantiles_to_plot=None, colors_quantiles=None, timewidth=0.5):
    """Plot effective population size (Ne) with confidence intervals and the other predictors"""

    # get the other predictors
    xml_file_path = os.path.join(path, 'xmls/HPAI_HLHxNx.glm.rep0.xml')

    with open(xml_file_path, 'r') as xml_file:
        xml_content = xml_file.read()

    # Extract actual predictors from the XML (only 3 exist: lpai, hpai, total)
    import re
    lpai_values = re.findall(r'<stateNode id="lpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    hpai_values = re.findall(r'<stateNode id="hpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    total_values = re.findall(r'<stateNode id="total" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)

    lpai_values = np.array([float(x) for x in lpai_values[0].split()])
    hpai_values = np.array([float(x) for x in hpai_values[0].split()])
    total_values = np.array([float(x) for x in total_values[0].split()])

    # set all values 0 to NaN
    lpai_values[lpai_values == 0] = np.nan
    hpai_values[hpai_values == 0] = np.nan
    total_values[total_values == 0] = np.nan

    # standardize the values
    lpai_values = (lpai_values - np.nanmean(lpai_values)) / np.nanstd(lpai_values)
    hpai_values = (hpai_values - np.nanmean(hpai_values)) / np.nanstd(hpai_values)
    total_values = (total_values - np.nanmean(total_values)) / np.nanstd(total_values)

    # Ensure arrays have the same length as valid_times
    target_length = len(valid_times)
    lpai_values = lpai_values[:target_length]
    hpai_values = hpai_values[:target_length]
    total_values = total_values[:target_length]

    # Plot confidence intervals as shaded ribbons for Ne
    if all_quantiles_ne is not None and quantiles_to_plot is not None and colors_quantiles is not None:
        # Define confidence intervals and their corresponding quantile pairs
        confidence_intervals = [
            (0.95, 0.025, 0.975),  # 95% CI: 2.5% to 97.5%
            (0.90, 0.05, 0.95),    # 90% CI: 5% to 95%
            (0.80, 0.1, 0.9),      # 80% CI: 10% to 90%
            (0.60, 0.2, 0.8),      # 60% CI: 20% to 80%
            (0.40, 0.3, 0.7),      # 40% CI: 30% to 70%
            (0.20, 0.4, 0.6),      # 20% CI: 40% to 60%
        ]
        
        # Plot each confidence interval as a shaded ribbon
        for ci_level, lower_q, upper_q in confidence_intervals:
            # Find the indices for these quantiles
            lower_idx = quantiles_to_plot.index(lower_q)
            upper_idx = quantiles_to_plot.index(upper_q)
            
            # Get the quantile values
            lower_values = all_quantiles_ne[:, lower_idx]
            upper_values = all_quantiles_ne[:, upper_idx]
            
            # Calculate alpha based on confidence level
            # Higher confidence = lighter shading, lower confidence = darker shading
            alpha = 0.1 + (1.0 - ci_level) * 0.4  # Range from 0.1 to 0.5
            
            # Create label for key confidence intervals
            label = f'{int(ci_level*100)}% CI' if ci_level in [0.95, 0.80, 0.40] else None
            
            # Create shaded ribbon for this confidence interval
            ax_ne.fill_between(valid_times, lower_values, upper_values,
                              alpha=alpha, color='#4A8B6F',  # Strong sage green
                              label=label, zorder=2, linewidth=0)

        # Plot the median as a thick shaded ribbon
        median_alpha = 0.6
        median_thickness = 0.01  # Small thickness for the median ribbon

        # Find median index (50% quantile)
        median_idx = quantiles_to_plot.index(0.5)
        median_values = all_quantiles_ne[:, median_idx]

        # Create upper and lower bounds for median ribbon
        median_upper = median_values + median_thickness
        median_lower = median_values - median_thickness

        ax_ne.fill_between(valid_times, median_lower, median_upper,
                          alpha=median_alpha, color='#4A8B6F',  # Strong sage green
                          label='Median', zorder=3, linewidth=0)
    else:
        # Fallback to original plotting if new parameters not provided
        ax_ne.fill_between(valid_times, lower_95_ne, upper_95_ne, alpha=0.2, color='#4E79A7', label='95% CI')
        ax_ne.fill_between(valid_times, lower_50_ne, upper_50_ne, alpha=0.4, color='#4E79A7', label='50% CI')
        ax_ne.plot(valid_times, medians_ne, color='#4E79A7', linewidth=2, label='Median')
    
    # Calculate the median Ne values for proper scaling
    if all_quantiles_ne is not None:
        median_idx = quantiles_to_plot.index(0.5)
        median_ne = all_quantiles_ne[:, median_idx]
    else:
        median_ne = medians_ne

    # Fit predictors to the Ne median using linear regression
    # This finds optimal scale and offset to match predictor variations to Ne variations
    from scipy import stats

    # Use LPAI as reference to determine scaling (both predictors will use same scale)
    valid_mask = ~np.isnan(lpai_values) & ~np.isnan(median_ne)
    # if np.sum(valid_mask) > 0:
    #     slope, intercept, _, _, _ = stats.linregress(lpai_values[valid_mask], median_ne[valid_mask])
    #     ne_scale = slope
    #     ne_offset = intercept
    # else:
    #     # Fallback if regression fails
    ne_range = np.max(median_ne) - np.min(median_ne)
    ne_scale = ne_range * 0.3
    ne_offset = np.mean(median_ne)

    ax_ne.plot(valid_times, lpai_values * ne_scale + ne_offset, color='#2ca02c', linewidth=1.5, label='LPAI', linestyle='--')
    ax_ne.plot(valid_times, hpai_values * ne_scale + ne_offset, color='#d62728', linewidth=1.5, label='HPAI', linestyle='--')

    # Add timeline shading to match the tree
    for i in np.arange(fromval, mrsi, 2 * timewidth):
        ax_ne.axvspan(i, i + timewidth, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)

    ax_ne.set_xlim(fromval, mrsi)

    # Set x-axis ticks to only show full years
    year_ticks = np.arange(np.ceil(fromval), np.floor(mrsi) + 1, 1.0)
    ax_ne.set_xticks(year_ticks)
    ax_ne.set_xticklabels([int(year) for year in year_ticks])

    # Remove top and right spines
    ax_ne.spines['top'].set_visible(False)
    ax_ne.spines['right'].set_visible(False)

    # Use real values on y-axis instead of log scale
    # Convert log Ne values back to real scale for display
    ax_ne.set_ylabel('Ne')

    # Set y-ticks with real values
    current_yticks = ax_ne.get_yticks()
    real_yticks = np.exp(current_yticks)
    ax_ne.set_yticklabels([f'{int(val):,}' if val >= 1 else f'{val:.2f}' for val in real_yticks])
  
def add_rate_axis(ax, time_points, offset_y, stretch_y, fromval, toval):
    """Add axis for reassortment rates with proper scaling"""
    axis_time = time_points[0]
    rate_range = np.log([0.0675, 0.125, 0.25, 0.5, 1, 2])
    axis_y_positions = rate_range * stretch_y + offset_y
    
    ax.plot([axis_time, axis_time], [axis_y_positions[0], axis_y_positions[-1]], 
            color='black', linewidth=1.5, zorder=10)
    
    for i, (rate_val, y_pos) in enumerate(zip(rate_range, axis_y_positions)):
        tick_length = (toval - fromval) * 0.005
        ax.plot([axis_time, axis_time + tick_length], [y_pos, y_pos], 
                color='black', linewidth=1, zorder=10)
        
        ax.plot([fromval, axis_time], [y_pos, y_pos], 
                color='grey', linewidth=0.5, alpha=1, linestyle='-', zorder=1)
        
        ax.text(axis_time + tick_length * 2, y_pos, f'{np.exp(rate_val):.3f}', 
                ha='left', va='center', fontsize=10, zorder=10)
    
    label_y = np.mean(axis_y_positions)
    ax.text(axis_time + (toval - fromval) * 0.05, label_y, 'Reassortment\nRate', 
            ha='center', va='center', fontsize=12, rotation=90, zorder=10)

def create_predictor_inset(ax, log_data):
    """Create bar plot for predictor support - now horizontal (flipped)"""
    ax_bar_inset = ax

    predictor_counts = log_data['predictorActive'].value_counts().sort_index()

    predictor_labels = {
        0: 'LPAI',
        1: 'HPAI', 
        2: 'Total AIV',
        3: 'Ne',
        4: 'Neither'
    }

    categories = list(range(len(predictor_labels)))
    counts = [predictor_counts.get(i, 0) for i in categories]
    labels = [predictor_labels[i] for i in categories]

    # Colors: LPAI, HPAI, Total AIV, Ne, Neither
    # Keep LPAI and HPAI colors to match other subplots, make others black/white
    colors = ['#2ca02c', '#d62728', '#000000', '#000000', '#000000']

    # Changed to bar() instead of barh() to make it vertical (flipped)
    bars = ax_bar_inset.bar(labels, [c/sum(counts) for c in counts], color=colors, alpha=0.8,
                            edgecolor='black', linewidth=0.2, width=0.7)

    ax_bar_inset.set_ylabel('Posterior Support for\nReassortment Rate Predictors')

    # Rotate x-axis labels for readability
    ax_bar_inset.tick_params(axis='x', rotation=45, labelsize=9)

    ax_bar_inset.spines['top'].set_visible(False)
    ax_bar_inset.spines['right'].set_visible(False)
    ax_bar_inset.spines['left'].set_linewidth(0.5)
    ax_bar_inset.spines['bottom'].set_linewidth(0.5)
    ax_bar_inset.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='y')
    ax_bar_inset.set_axisbelow(True)
    ax_bar_inset.set_ylim(0, 0.5)

    ax_bar_inset.patch.set_facecolor('white')
    ax_bar_inset.patch.set_alpha(0.9)

    for spine in ax_bar_inset.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(0.5)

def create_clade_probability_plot(ax_clade_probs, ax_clade_heights, mrsi, path, fromval, toval):
    """Create clade probability and height plots"""
    log_file_path = os.path.join(path, 'combined/HPAI_HLHxNx.glm.cladeprobs.csv')
    log_data = pd.read_csv(log_file_path, sep=',')

    # Create boxplot for reassortment event probabilities by clade type
    clade_probs = log_data['no_event_prob'].values
    clade_types = log_data['clade'].values
    clade_types_unique = np.unique(clade_types)
    clade_probs_by_type = {clade: [] for clade in clade_types_unique}

    for prob, clade in zip(clade_probs, clade_types):
        clade_probs_by_type[clade].append(prob)

    # Create boxplot for probabilities without outliers
    ax_clade_probs.boxplot(clade_probs_by_type.values(), labels=clade_probs_by_type.keys(), showfliers=False)
    ax_clade_probs.set_ylabel('Probability of at least one\nreassortment event')
    ax_clade_probs.tick_params(axis='x', rotation=0)

    # Style the probability plot
    ax_clade_probs.spines['top'].set_visible(False)
    ax_clade_probs.spines['right'].set_visible(False)
    ax_clade_probs.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='y')

    # Only create clade heights plot if ax_clade_heights is provided
    if ax_clade_heights is not None:
        # Calculate clade heights (time from MRSI to minimum time in clade)
        clade_heights = mrsi - log_data['min_time'].values
        clade_heights_by_type = {clade: [] for clade in clade_types_unique}

        for height, clade in zip(clade_heights, clade_types):
            clade_heights_by_type[clade].append(height)

        # Add timeline shading to match the other plots
        for i in np.arange(fromval, toval, 2 * 1):
            ax_clade_heights.axvspan(i, i + 1, facecolor='#F2F2F2', edgecolor='none', alpha=1, zorder=0)

        # Create horizontal violin plot for clade heights (violins along x-axis)
        clade_names = list(clade_heights_by_type.keys())
        y_positions = np.arange(len(clade_names))

        parts = ax_clade_heights.violinplot(
            clade_heights_by_type.values(),
            positions=y_positions,
            showmeans=False,
            showmedians=True,
            widths=0.8,
            vert=False  # This makes violins horizontal (along x-axis)
        )

        # Style the violin plot
        for pc in parts['bodies']:
            pc.set_facecolor('#4E79A7')
            pc.set_alpha(0.6)
            pc.set_edgecolor('black')
            pc.set_linewidth(0.5)

        # Style median lines
        if 'cmedians' in parts:
            parts['cmedians'].set_color('black')
            parts['cmedians'].set_linewidth(1.5)

        ax_clade_heights.set_xlabel('Clade Height (years)')
        ax_clade_heights.set_yticks(y_positions)
        ax_clade_heights.set_yticklabels(clade_names)

        # Set x-axis limits to match the timeline
        ax_clade_heights.set_xlim(fromval, toval)

        # Style the heights plot
        ax_clade_heights.spines['top'].set_visible(False)
        ax_clade_heights.spines['right'].set_visible(False)
        ax_clade_heights.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='x')

def create_reassortment_probability_plot(ax, path):
    """Create reassortment probability plot from normalized co-infection data"""
    log_file_path = os.path.join(path, 'SIR_SIS/co_inf_normalized.txt')
    log_data = pd.read_csv(log_file_path, sep='\t')

    # Extract data
    time = log_data['time'].values
    counts = log_data['normalized_counts'].values
    method = log_data['method'].values
    
    # Get unique methods
    unique_methods = np.unique(method)
    
    # Define colors for the two methods
    method_colors = {
        'transmission rate limited': '#d62728',  # Red
        'susceptible depletion': '#2ca02c'       # Green
    }
    
    # Plot line for each method
    for method_name in unique_methods:
        mask = method == method_name
        method_time = time[mask]
        method_counts = counts[mask]
        
        # Sort by time to ensure proper line plotting
        sort_idx = np.argsort(method_time)
        method_time_sorted = method_time[sort_idx]
        method_counts_sorted = method_counts[sort_idx]
        
        ax.plot(method_time_sorted, method_counts_sorted, 
                color=method_colors.get(method_name, '#4E79A7'),
                linewidth=2, 
                label=method_name)
    
    
    # Set axis properties
    ax.set_xlim(-3, 3)
    ax.set_xlabel('Time relative to peak')
    ax.set_ylabel('Normalized Co-infection Counts')
    ax.legend(fontsize=8)
    
    # Style the plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.3)

def create_cluster_size_comparison_plot(ax, path):
    """Create plot showing probability that children with reassortment > children without, across posterior"""
    print("DEBUG: Starting create_cluster_size_comparison_plot")
    cluster_file_path = os.path.join(path, 'combined/HPAI_HLHxNx.glm.cluster_comparison.txt')
    print(f"DEBUG: Looking for file at {cluster_file_path}")
    
    try:
        cluster_data = pd.read_csv(cluster_file_path, sep='\t')
        
        # Group by iteration and calculate probability for each iteration
        iteration_probs = []
        
        for iteration in cluster_data['iteration'].unique():
            iter_data = cluster_data[cluster_data['iteration'] == iteration]
            
            # Calculate probability that leafsWith > leafsWithout for this iteration
            with_greater = (iter_data['leafsWith'] > iter_data['leafsWithout']).sum()
            total_comparisons = len(iter_data)
            prob_with_greater = with_greater / total_comparisons
            
            iteration_probs.append(prob_with_greater)

        # Plot histogram of probabilities across posterior
        n, bins, patches = ax.hist(iteration_probs, bins=30, alpha=0.7, color='#1f77b4',
                density=True, edgecolor='black', linewidth=0.5)

        ax.set_xlabel('Probability (With > Without)')
        ax.set_ylabel('Density')
        # Title removed as requested

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='y')
        ax.set_xlim(0, 1)

        # Debug: print some info
        print(f"Cluster plot: {len(iteration_probs)} iterations, max_bin_height={max(n):.3f}")
        
    except FileNotFoundError:
        ax.text(0.5, 0.5, 'Cluster comparison file not found', 
                transform=ax.transAxes, ha='center', va='center')
        ax.set_title('Posterior Distribution of\nP(Children With > Children Without)')
    except Exception as e:
        ax.text(0.5, 0.5, f'Error loading cluster data:\n{str(e)}', 
                transform=ax.transAxes, ha='center', va='center', fontsize=8)
        ax.set_title('Posterior Distribution of\nP(Children With > Children Without)')

def save_figure(fig, approach):
    """Save the figure with appropriate filename"""
    # save the figure into the /Users/nmueller/Documents/github/CoInfection-Material/Figures directory
    figname = 'H5N1_' + approach + '_combined.pdf'
    figname = os.path.join('/Users/nmueller/Documents/github/CoInfection-Material/Figure1.pdf')
    print('Saving figure to:')
    print(figname)
    fig.savefig(figname, bbox_inches='tight', pad_inches=0.2)
    

def main(force=False):
    """Main function to run the complete analysis"""
    # Setup
    setup_matplotlib()
    colours, segments, new_order, colour_cycle, path, linewidth, approach = define_constants()
    
    # Process GLM data (replaces convert_GLM_data.R)
    process_glm_data(path, force=force)
    
    # Create figure with better spacing
    fig = plt.figure(figsize=(18, 12), facecolor='w')

    # Use 20 columns for finer control and equal widths
    # Set different wspace for different column pairs to add more padding between E and F
    gs = fig.add_gridspec(3, 20,
                          width_ratios=[1]*20,
                          height_ratios=[1, 1, 1],
                          hspace=0.35)  # Vertical spacing between rows

    mrsi = 2025.12877
    start_time = 2021.5

    # Network takes rows 0-1, columns 0-11 (11/20 = 55%)
    ax = fig.add_subplot(gs[0:2, :11], facecolor='w')
    # Predictor plot - top right (columns 13-19, 7 columns wide) with gap at column 12
    ax_probs = fig.add_subplot(gs[0, 13:], facecolor='w')
    # Two smaller plots in middle right - equally wide (3 columns each, with gaps)
    ax_clade_probs = fig.add_subplot(gs[1, 13:16], facecolor='w')  # 3 columns
    ax_cluster_comparison = fig.add_subplot(gs[1, 17:], facecolor='w')  # 3 columns (16 is gap)
    # Reassortment and Ne plots - with small gap between them, filling entire bottom row
    ax_rates = fig.add_subplot(gs[2, :9], facecolor='w')  # columns 0-8
    ax_ne = fig.add_subplot(gs[2, 10:], facecolor='w')  # columns 10-19 (column 9 is gap)

    # Process tree
    ll = load_and_process_tree(bt, path, approach)
    fromval, toval, timewidth = setup_tree_plot(ax, ll)

    # Initialize and process tree data
    initialize_traits(ll)
    assign_reassortment_colors(ll, colour_cycle)
    reassortment_events = draw_tree_branches(ax, ll, colour_cycle, linewidth, fromval)
    finalize_tree_plot(ax, ll, fromval, toval)

    # Add subplot labels - positioned outside and above each subplot
    # A and B should align at the same vertical position (top of row 0)
    ax.text(-0.02, 1.02, 'A', transform=ax.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_probs.text(-0.02, 1.02, 'B', transform=ax_probs.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_clade_probs.text(-0.02, 1.02, 'C', transform=ax_clade_probs.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_cluster_comparison.text(-0.02, 1.02, 'D', transform=ax_cluster_comparison.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_rates.text(-0.02, 1.02, 'E', transform=ax_rates.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')
    ax_ne.text(-0.02, 1.02, 'F', transform=ax_ne.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')
    
    # Process log data and plot rates
    log_data = load_log_data(path, approach)
    (valid_times, medians, lower_95, upper_95, lower_50, upper_50,
     medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, 
     time_points, all_quantiles, all_quantiles_ne, quantiles_to_plot, colors_quantiles) = calculate_rate_quantiles(log_data, mrsi, path)
    offset_y, stretch_y = plot_reassortment_rates(ax_rates, valid_times, medians, lower_95, upper_95, lower_50, upper_50, start_time, mrsi, path, all_quantiles, quantiles_to_plot, colors_quantiles, timewidth)
    plot_ne(ax_ne, valid_times, medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, start_time, mrsi, path, all_quantiles_ne, quantiles_to_plot, colors_quantiles, timewidth)
    
    # Create inset and finalize
    create_predictor_inset(ax_probs, log_data)
    create_clade_probability_plot(ax_clade_probs, None, mrsi, path, 2023, toval)
    create_cluster_size_comparison_plot(ax_cluster_comparison, path)
    
    print(toval)
    # Tree is now rotated 90 degrees, so don't override the xlim
    # ax.set_xlim(fromval, toval)

    # Adjust subplot positions to prevent y-axis label overlap
    # Using subplots_adjust for precise control over spacing
    fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.08, hspace=0.35, wspace=0.25)
    
    save_figure(fig, approach)
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate H5N1 GLM analysis plots')
    parser.add_argument('--force', action='store_true', 
                       help='Force rerun of all BEAST commands and calculations, even if output files exist')
    
    args = parser.parse_args()
    main(force=args.force)