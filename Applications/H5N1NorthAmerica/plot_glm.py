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
    colour_cycle = ['#969696', '#737373', '#525252', '#252525']
    path = '/Users/nmueller/Documents/github/CoInfection-Material/Applications/H5N1NorthAmerica/'
    linewidth = 0.5
    approach = 'glm'
    
    return colours, segments, new_order, colour_cycle, path, linewidth, approach


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
    [ax.spines[loc].set_visible(False) for loc in ax.spines if loc != 'bottom']
    
    fromval = float(int(ll.root.absoluteTime + ll.treeHeight)) - 5
    toval = float(int(ll.root.absoluteTime + ll.treeHeight)) + 0.5
    timewidth = 0.5
    
    for i in np.arange(fromval, toval, 2 * timewidth):
        ax.axvspan(i, i + timewidth, facecolor='#F2F2F2', edgecolor='none', alpha=1, zorder=0)
    
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
    """Draw all tree branches and nodes"""
    reassortment_events = []
    
    for k in ll.Objects:
        x = k.absoluteTime
        xp = k.parent.absoluteTime
        if xp != None:
            xp = max(xp, fromval + 0.000001)
        y = k.y
        col = colour_cycle[k.traits['re'] % len(colour_cycle)]
        col = 'black'
        
        if isinstance(k, bt.reticulation) == False:
            col_lin = 'grey'
            lw_scale = 1
            if k.traits['seg0'] == 'true':
                col_lin = 'black'
                lw_scale = 2
            
            ax.plot([x, xp], [y, y], color=col_lin, lw=linewidth * lw_scale, solid_capstyle='round')
        else:
            ax.plot([x, xp], [y, y], color=col, lw=linewidth, ls=':', solid_capstyle='round')
        
        if k.branchType == 'node':
            left, right = k.children[-1].y, k.children[0].y
            
            col_lin1 = 'grey'
            lw_scale1 = 1
            col_lin2 = 'grey'
            lw_scale2 = 1
            
            if k.children[-1].traits['seg0'] == 'true':
                col_lin1 = 'black'
                lw_scale1 = 2
            
            if k.children[0].traits['seg0'] == 'true':
                col_lin2 = 'black'
                lw_scale2 = 2
            
            ax.plot([x, x], [left, k.y], color=col_lin1, lw=linewidth * lw_scale1, solid_capstyle='round')
            ax.plot([x, x], [k.y, right], color=col_lin2, lw=linewidth * lw_scale2, solid_capstyle='round')
        
        elif isinstance(k, bt.leaf):
            ax.scatter(x, y, s=30, facecolor='black', edgecolor='none', zorder=4)
            ax.scatter(x, y, s=10, facecolor='grey', edgecolor='none', zorder=5)
        
        elif isinstance(k, bt.reticulation):
            segs = sorted(map(int, k.traits['segments']))
            
            reassortment_events.append({
                'time': x,
                'segments': segs,
                'posterior': k.traits.get('posterior', 1.0)
            })
            
            ax.scatter(x, k.target.y, s=20, facecolor="black", edgecolor='none', zorder=4)
            ax.scatter(x, k.target.y, s=10, facecolor="#d62728", edgecolor='none', zorder=5)
            ax.plot([x, x], [y, k.target.y], color=col, lw=linewidth, ls='-', solid_capstyle='round')
            
            for i in range(len(segs)):
                name = segs[i]
                c = 'black'
                o = 1 / 20.
                posterior_val = round(k.traits['posterior'], 2)
    
    return reassortment_events

def finalize_tree_plot(ax, ll, fromval, toval):
    """Apply final formatting to the tree plot"""
    ax.set_yticks([])
    ax.set_ylim(ll.ySpan * 1.01, -ll.ySpan * 0.05)
    ax.set_xlim(fromval, toval)
    ax.spines['bottom'].set_visible(False)
    # ax.set_xlabel('Time')

def load_log_data(path, approach):
    """Load and process the log file data"""
    log_file_path = os.path.join(path, 'combined/HPAI_HLHxNx.' + approach + '.log')
    print(sys.executable)
    print(log_file_path)
    log_data = pd.read_csv(log_file_path, sep='\t')
    return log_data

def calculate_rate_quantiles(log_data, mrsi):
    """Calculate quantiles for reassortment rates over time"""
    rate_shifts = np.linspace(0, 4.13150684931507, 40)
    
    time_points = [mrsi - shift for shift in rate_shifts]
    
    quantiles_to_plot = [0.05, 0.25, 0.5, 0.75, 0.95]
    colors_quantiles = ['#E8F4FD', '#B3D9F2', '#4E79A7', '#B3D9F2', '#E8F4FD']
    
    medians = []
    lower_95 = []
    upper_95 = []
    lower_50 = []
    upper_50 = []
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
        
        medians.append(np.quantile(rates_at_time, 0.5))
        lower_95.append(np.quantile(rates_at_time, 0.025))
        upper_95.append(np.quantile(rates_at_time, 0.975))
        lower_50.append(np.quantile(rates_at_time, 0.25))
        upper_50.append(np.quantile(rates_at_time, 0.75))
        valid_times.append(time_points[i])
        
        medians_ne.append(np.quantile(ne_at_time, 0.5))
        lower_95_ne.append(np.quantile(ne_at_time, 0.025))
        upper_95_ne.append(np.quantile(ne_at_time, 0.975))
        lower_50_ne.append(np.quantile(ne_at_time, 0.25))
        upper_50_ne.append(np.quantile(ne_at_time, 0.75))
    
    return (valid_times, medians, lower_95, upper_95, lower_50, upper_50,
            medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, time_points)

def plot_reassortment_rates(ax, valid_times, medians, lower_95, upper_95, lower_50, upper_50, fromval, toval, path):
    """Plot reassortment rates with confidence intervals"""
    offset_y = 500.
    stretch_y = -25
    xml_file_path = os.path.join(path, 'xmls/HPAI_HLHxNx.glm.rep0.xml')

    with open(xml_file_path, 'r') as xml_file:
        xml_content = xml_file.read()

    # Extract the relevant lines for the predictors
    import re

    lpai_values = re.findall(r'<stateNode id="lpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    hpai_values = re.findall(r'<stateNode id="hpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    total_values = re.findall(r'<stateNode id="total" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    overlap_values = re.findall(r'<stateNode id="overlap" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)   

    lpai_values = np.array([float(x) for x in lpai_values[0].split()])
    hpai_values = np.array([float(x) for x in hpai_values[0].split()])
    total_values = np.array([float(x) for x in total_values[0].split()])
    overlap_values = np.array([float(x) for x in overlap_values[0].split()])

    # set all values 0 to NaN
    lpai_values[lpai_values == 0] = np.nan
    hpai_values[hpai_values == 0] = np.nan
    total_values[total_values == 0] = np.nan
    overlap_values[overlap_values == 0] = np.nan

    # standardize the values
    lpai_values = (lpai_values - np.nanmean(lpai_values)) / np.nanstd(lpai_values)
    hpai_values = (hpai_values - np.nanmean(hpai_values)) / np.nanstd(hpai_values)
    total_values = (total_values - np.nanmean(total_values)) / np.nanstd(total_values)
    overlap_values = (overlap_values - np.nanmean(overlap_values)) / np.nanstd(overlap_values)

    
    medians = np.array(medians)
    lower_95 = np.array(lower_95)
    upper_95 = np.array(upper_95)
    lower_50 = np.array(lower_50)
    upper_50 = np.array(upper_50)
    
    ax.fill_between(valid_times, lower_95, upper_95, alpha=0.2, color='#4E79A7', label='95% CI', zorder=2)
    ax.fill_between(valid_times, lower_50, upper_50, alpha=0.4, color='#4E79A7', label='50% CI', zorder=2)

    ax.plot(valid_times, medians, color='#4E79A7', linewidth=2, label='Median', zorder=2)

    # Plot the other predictors
    ax.plot(valid_times, lpai_values, color='#d62728', linewidth=1.5, label='LPAI', linestyle='--')
    ax.plot(valid_times, hpai_values, color='#ff7f0e', linewidth=1.5, label='HPAI', linestyle='--')

    # Add timeline shading to match the tree
    for i in np.arange(fromval, toval, 2 * 1):
        ax.axvspan(i, i + 1, facecolor='#F2F2F2', edgecolor='none', alpha=1, zorder=0)

    ax.set_xlim(fromval, toval)
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


def plot_ne(ax_ne, valid_times, medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, fromval, toval, path):
    """Plot effective population size (Ne) with confidence intervals and the other predictors"""

    # get the other predictors
    xml_file_path = os.path.join(path, 'xmls/HPAI_HLHxNx.glm.rep0.xml')

    with open(xml_file_path, 'r') as xml_file:
        xml_content = xml_file.read()

    # Extract the relevant lines for the predictors
    import re
    lpai_values = re.findall(r'<stateNode id="lpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    hpai_values = re.findall(r'<stateNode id="hpai" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    total_values = re.findall(r'<stateNode id="total" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)
    overlap_values = re.findall(r'<stateNode id="overlap" spec="parameter.RealParameter" value="([^"]+)"/>', xml_content)   
    lpai_values = np.array([float(x) for x in lpai_values[0].split()])
    hpai_values = np.array([float(x) for x in hpai_values[0].split()])
    total_values = np.array([float(x) for x in total_values[0].split()])
    overlap_values = np.array([float(x) for x in overlap_values[0].split()])

    # set all values 0 to NaN
    lpai_values[lpai_values == 0] = np.nan
    hpai_values[hpai_values == 0] = np.nan
    total_values[total_values == 0] = np.nan
    overlap_values[overlap_values == 0] = np.nan

    # standardize the values
    lpai_values = (lpai_values - np.nanmean(lpai_values)) / np.nanstd(lpai_values)
    hpai_values = (hpai_values - np.nanmean(hpai_values)) / np.nanstd(hpai_values)
    total_values = (total_values - np.nanmean(total_values)) / np.nanstd(total_values)
    overlap_values = (overlap_values - np.nanmean(overlap_values)) / np.nanstd(overlap_values)

    ax_ne.fill_between(valid_times, lower_95_ne, upper_95_ne, alpha=0.2, color='#4E79A7', label='95% CI')
    ax_ne.fill_between(valid_times, lower_50_ne, upper_50_ne, alpha=0.4, color='#4E79A7', label='50% CI')
    ax_ne.plot(valid_times, medians_ne, color='#4E79A7', linewidth=2, label='Median')
    
    print(overlap_values)
    print(valid_times)
    # Plot the other predictors
    ax_ne.plot(valid_times, lpai_values+3, color='#d62728', linewidth=1.5, label='LPAI', linestyle='--')
    ax_ne.plot(valid_times, hpai_values+3, color='#ff7f0e', linewidth=1.5, label='HPAI', linestyle='--')
    # ax_ne.plot(valid_times, total_values+3, color='#2ca02c', linewidth=1.5, label='Total AIV', linestyle='--')
    # ax_ne.plot(valid_times, overlap_values+3, color='#1f77b4', linewidth=1.5, label='Overlap', linestyle='--')

    # Add timeline shading to match the tree
    for i in np.arange(fromval, toval, 2 * 1):
        ax_ne.axvspan(i, i + 1, facecolor='#F2F2F2', edgecolor='none', alpha=1, zorder=0)

    ax_ne.set_xlim(fromval, toval)
    # ax_ne.set_xlabel('Time')
    ax_ne.set_ylabel('Ne')

    
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
    """Create inset bar plot for predictor support"""
    ax_bar_inset = ax
    
    predictor_counts = log_data['predictorActive'].value_counts().sort_index()
    
    predictor_labels = {
        0: 'H5 LPAI',
        1: 'LPAI',
        2: 'HPAI', 
        3: 'Total AIV',
        4: 'Overlap',
        5: 'Ne',
        6: 'Neither'
    }
    
    categories = list(range(len(predictor_labels)))
    counts = [predictor_counts.get(i, 0) for i in categories]
    labels = [predictor_labels[i] for i in categories]
    
    colors = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4', '#9467bd', '#8c564b']
    bars = ax_bar_inset.barh(labels, counts/sum(counts), color=colors, alpha=0.8, 
                            edgecolor='black', linewidth=0.2, height=0.7)
    
    ax_bar_inset.set_xlabel('Posterior Support for\nReassortment Rate Predictors')
    
    ax_bar_inset.xaxis.set_ticks_position('top')
    ax_bar_inset.xaxis.set_label_position('top')
    
    ax_bar_inset.spines['top'].set_visible(False)
    ax_bar_inset.spines['right'].set_visible(False)
    ax_bar_inset.spines['left'].set_linewidth(0.5)
    ax_bar_inset.spines['bottom'].set_linewidth(0.5)
    ax_bar_inset.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='x')
    ax_bar_inset.set_axisbelow(True)
    ax_bar_inset.set_xlim(0, 1)
    
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

    # Create boxplot for probabilities
    ax_clade_probs.boxplot(clade_probs_by_type.values(), labels=clade_probs_by_type.keys())
    ax_clade_probs.set_ylabel('Probability of at least one\nreassortment event')
    ax_clade_probs.tick_params(axis='x', rotation=0)
    
    # Style the probability plot
    ax_clade_probs.spines['top'].set_visible(False)
    ax_clade_probs.spines['right'].set_visible(False)
    ax_clade_probs.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='y')

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


def save_figure(fig, approach):
    """Save the figure with appropriate filename"""
    # save the figure into the /Users/nmueller/Documents/github/CoInfection-Material/Figures directory
    figname = 'H5N1_' + approach + '_combined.pdf'
    figname = os.path.join('/Users/nmueller/Documents/github/CoInfection-Material/Figure1.pdf')
    print('Saving figure to:')
    print(figname)
    fig.savefig(figname, bbox_inches='tight', pad_inches=0.2)
    

def main():
    """Main function to run the complete analysis"""
    # Setup
    setup_matplotlib()
    colours, segments, new_order, colour_cycle, path, linewidth, approach = define_constants()
    
    # Create figure with better spacing
    fig = plt.figure(figsize=(16, 12), facecolor='w')
    
    gs = fig.add_gridspec(3, 4,
                          width_ratios=[5, 1, 1, 1],
                          height_ratios=[1, 1, 1],
                          hspace=0.4,  # Increased vertical spacing
                          wspace=0.3)  # Increased horizontal spacing
    
    mrsi = 2025.12877
    start_time = 2021
    
    ax = fig.add_subplot(gs[:, 0], facecolor='w')
    ax_clade_probs = fig.add_subplot(gs[2, 2], facecolor='w')
    ax_clade_heights = fig.add_subplot(gs[2, 3], facecolor='w')
    ax_probs = fig.add_subplot(gs[2, 1], facecolor='w')
    ax_ne = fig.add_subplot(gs[1, 1:], facecolor='w')
    ax_rates = fig.add_subplot(gs[0, 1:], facecolor='w')

    # Process tree
    ll = load_and_process_tree(bt, path, approach)
    fromval, toval, timewidth = setup_tree_plot(ax, ll)
    
    # Initialize and process tree data
    initialize_traits(ll)
    assign_reassortment_colors(ll, colour_cycle)
    reassortment_events = draw_tree_branches(ax, ll, colour_cycle, linewidth, fromval)
    finalize_tree_plot(ax, ll, fromval, toval)
    
    # Process log data and plot rates
    log_data = load_log_data(path, approach)
    (valid_times, medians, lower_95, upper_95, lower_50, upper_50,
     medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, time_points) = calculate_rate_quantiles(log_data, mrsi)
    offset_y, stretch_y = plot_reassortment_rates(ax_rates, valid_times, medians, lower_95, upper_95, lower_50, upper_50, start_time, toval, path)
    plot_ne(ax_ne, valid_times, medians_ne, lower_95_ne, upper_95_ne, lower_50_ne, upper_50_ne, start_time, toval, path)
    
    # Create inset and finalize
    create_predictor_inset(ax_probs, log_data)
    create_clade_probability_plot(ax_clade_probs, ax_clade_heights, mrsi, path, 2023, toval)
    # reset everthing about ax_clade_heights
    ax_clade_heights.clear()
    create_reassortment_probability_plot(ax_clade_heights, path)
    
    print(toval)
    ax.set_xlim(fromval, toval)
    
    # Use constrained_layout for better automatic spacing
    plt.tight_layout(pad=2.0)  # Add extra padding
    
    save_figure(fig, approach)
    plt.show()
    plt.close(fig)

if __name__ == "__main__":
    main()