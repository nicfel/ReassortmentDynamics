import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import os
import sys
import random
from datetime import datetime, timedelta
from matplotlib.patches import FancyBboxPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Add baltic_bacter path
sys.path.append('/Users/nmueller/Documents/github/CoInfection-Material/Applications/NetworkViz/')
import baltic as bt


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
    rate_shift_str = '0 1.25 2.5 3.75 5 6 9.5 13 16.5 20'
    rate_shifts = np.array([float(x) for x in rate_shift_str.split()])
    
    # Set random seed for reproducibility
    np.random.seed(6465546)
    
    # Define dates
    mrsi = pd.to_datetime("2025-02-17")
    mrsi_hpai = mrsi
    mrsi_lpai = pd.to_datetime("2024-08-22")
    
    # Convert to decimal year for baltic
    mrsi_dec = 2025 + (mrsi - pd.to_datetime("2025-01-01")).days / 365.25
    
    # Define segment order
    segment_order = ["HA", "NA", "MP", "NS", "NP", "PB1", "PB2", "PA"]
    
    # Define color schemes
    methods_colors = {
        "constant": "#E41A1C",      # red
        "skyline": "#377EB8",       # blue
        "skyline_Ne": "#4DAF4A"     # green
    }
    
    species_colors = {
        "cow": "#238B45",   # teal
        "bird": "#D95F02"   # burnt orange
    }
    
    clade_colors = {
        "B3.13": "#238B45",  # teal
        "D1.1": "#D95F02"    # burnt orange
    }
    
    lineage_colors = {
        "HPAI": "#E41A1C",    # red
        "LPAI": "#377EB8",    # blue
        "unknown": "#4DAF4A"  # green
    }
    
    return (clades, rate_shifts, mrsi, mrsi_hpai, mrsi_lpai, mrsi_dec,
            segment_order, methods_colors, species_colors, 
            clade_colors, lineage_colors)


def load_and_process_tree(path, segment, is_constant=True):
    """Load and process the phylogenetic tree using baltic_bacter"""
    if is_constant:
        tree_file = f'combined/HLHxNx.constant.{segment}.tree'
    else:
        tree_file = f'combined/HLHxNx.dependent.{segment}.tree'
    
    tree_path = os.path.join(path, tree_file)
    
    if not os.path.exists(tree_path):
        print(f"Warning: Tree file not found: {tree_path}")
        return None
    
    print(f"Loading tree: {tree_path}")
    ll = bt.loadNexus(tree_path, date_fmt='%Y-%m-%d', verbose=True)
    ll.drawTree()
    return ll


def setup_tree_axes(ax, ll, fromval=2020, toval=2025.5):
    """Configure the tree plot axes and timeline"""
    ax.set_facecolor('w')
    [ax.spines[loc].set_visible(False) for loc in ax.spines if loc != 'bottom']
    
    # Add alternating background shading
    timewidth = 0.5
    for i in np.arange(fromval, toval, 2 * timewidth):
        ax.axvspan(i, i + timewidth, facecolor='#F2F2F2', edgecolor='none', alpha=1, zorder=0)
    
    return fromval, toval, timewidth


def draw_ha_tree_with_clades(ax, ll, mrsi_dec, linewidth=0.5):
    """Draw HA tree with HPAI/LPAI clade coloring"""
    if ll is None:
        ax.text(0.5, 0.5, "Tree data not available", 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    fromval, toval, _ = setup_tree_axes(ax, ll, 2020, mrsi_dec + 0.05)
    
    # Draw branches colored by HPAI+LPAI trait
    for k in ll.Objects:
        x = k.absoluteTime
        if k.parent:
            xp = k.parent.absoluteTime
            print(xp)
            xp = max(xp, fromval + 0.000001)
        else:
            xp = x
            
        y = k.y
        
        # Get HPAI+LPAI trait value for coloring
        clade_support = 0
        if hasattr(k, 'traits') and 'HPAI+LPAI' in k.traits:
            try:
                clade_support = float(k.traits['HPAI+LPAI'])
            except (ValueError, TypeError):
                clade_support = 0
        
        # Color based on support value
        if clade_support > 0.7:
            color = '#F93946'
        else:
            color = 'black'
        
        # Draw horizontal branch
        if k.branchType == 'leaf':
            ax.plot([x, xp], [y, y], color=color, lw=linewidth, solid_capstyle='round')
            
            # Check if it's a cow tip
            is_cow = 'cow' in k.name.lower() if hasattr(k, 'name') else False
            if is_cow:
                ax.scatter(x, y, s=30, facecolor='black', edgecolor='none', zorder=4)
                ax.scatter(x, y, s=15, facecolor='#238B45', edgecolor='none', zorder=5)
            else:
                ax.scatter(x, y, s=20, facecolor='black', edgecolor='none', zorder=4)
                ax.scatter(x, y, s=10, facecolor='grey60', edgecolor='none', zorder=5)
        
        elif k.branchType == 'node':
            ax.plot([x, xp], [y, y], color=color, lw=linewidth, solid_capstyle='round')
            
            # Draw vertical lines to children
            if hasattr(k, 'children') and len(k.children) >= 2:
                left = k.children[-1].y
                right = k.children[0].y
                ax.plot([x, x], [left, right], color=color, lw=linewidth, solid_capstyle='round')
    
    # Format axes
    ax.set_xlim(fromval, toval)
    ax.set_ylim(ll.ySpan * 1.01, -ll.ySpan * 0.05)
    ax.set_xticks(range(2020, 2026))
    ax.set_xlabel('Year')
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)


def draw_segment_tree(ax, ll, mrsi_dec, lineage_colors, segment_name, linewidth=0.5):
    """Draw tree for non-HA segments with HPAI/LPAI tip coloring"""
    if ll is None:
        ax.text(0.5, 0.5, f"{segment_name} tree\nnot available", 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    fromval, toval, _ = setup_tree_axes(ax, ll, 2020, mrsi_dec + 0.05)
    
    # Draw all branches
    for k in ll.Objects:
        x = k.absoluteTime
        if k.parent:
            xp = k.parent.absoluteTime
            xp = max(xp, fromval + 0.000001)
        else:
            xp = x
        y = k.y
        
        # Default color
        branch_color = 'grey'
        
        # Draw horizontal branch
        if k.branchType == 'leaf':
            ax.plot([x, xp], [y, y], color=branch_color, lw=linewidth, solid_capstyle='round')
            
            # Determine tip type based on trait or name
            tip_type = 'unknown'
            if hasattr(k, 'traits') and 'type' in k.traits:
                tip_type = k.traits['type']
            elif hasattr(k, 'name'):
                # Simple heuristic: check if HPAI/LPAI in name
                if 'HPAI' in k.name or 'H5N1' in k.name:
                    tip_type = 'HPAI'
                elif 'LPAI' in k.name:
                    tip_type = 'LPAI'
            
            # Draw tip point
            tip_color = lineage_colors.get(tip_type, '#4DAF4A')
            ax.scatter(x, y, s=20, facecolor='black', edgecolor='none', zorder=4)
            ax.scatter(x, y, s=10, facecolor=tip_color, edgecolor='none', zorder=5)
        
        elif k.branchType == 'node':
            ax.plot([x, xp], [y, y], color=branch_color, lw=linewidth, solid_capstyle='round')
            
            # Draw vertical lines to children
            if hasattr(k, 'children') and len(k.children) >= 2:
                left = k.children[-1].y
                right = k.children[0].y
                ax.plot([x, x], [left, right], color=branch_color, lw=linewidth, solid_capstyle='round')
    
    # Format axes
    ax.set_xlim(fromval, toval)
    ax.set_ylim(ll.ySpan * 1.01, -ll.ySpan * 0.05)
    ax.set_xticks(range(2020, 2026, 2))  # Every 2 years
    ax.set_xlabel('')
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)


def load_log_file(path, is_constant=True):
    """Load and process log files"""
    if is_constant:
        log_file = pd.read_csv(f"{path}/combined/HLHxNx.constant.log", sep='\t')
    else:
        log_file = pd.read_csv(f"{path}/combined/HLHxNx.dependent.log", sep='\t')
    return log_file


def load_cases_data(path):
    """Load and process positive cases data"""
    cases = pd.read_csv(f"{path}/tables/APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv")
    cases['date'] = pd.to_datetime(cases['Date_Collected'], format='%Y-%m-%d')
    return cases


def smooth_case_data(cases):
    """Calculate smoothed case data with moving average"""
    min_date = cases['date'].min()
    max_date = cases['date'].max()
    smoothed_data = []
    
    for d in pd.date_range(min_date + pd.Timedelta(days=23), 
                           max_date - pd.Timedelta(days=23), 
                           freq='D'):
        # Get all instances within time window
        window = cases[(cases['date'] >= d - pd.Timedelta(days=23)) & 
                      (cases['date'] <= d + pd.Timedelta(days=23))]
        
        # Calculate total AIV detected
        total_AIV = (window['Final_IAV'] == "Detected").sum()
        
        # Calculate high path cases
        total_HPAI = ((window['Final_H5'] == "Detected") & 
                     (window['Final_Pathogenicity'] == "High Path AI")).sum()
        
        # LPAI data
        smoothed_data.append({
            'date': d,
            'positivity': (total_AIV - total_HPAI) / len(window) if len(window) > 0 else 0,
            'type': 'LPAI'
        })
        
        # HPAI data
        smoothed_data.append({
            'date': d,
            'positivity': total_HPAI / len(window) if len(window) > 0 else 0,
            'type': 'HPAI'
        })
    
    return pd.DataFrame(smoothed_data)


def calculate_reassortment_rates(log_file, rate_shifts, mrsi, is_constant=True):
    """Calculate reassortment rate quantiles over time"""
    reassortment_data = []
    
    for i, shift in enumerate(rate_shifts):
        col_name = f"InfectedToRho.{i}"
        if col_name not in log_file.columns:
            continue
            
        rate = log_file[col_name].values
        
        if not is_constant:
            ne_col = f"logNe.{i}"
            if ne_col in log_file.columns:
                rate = rate + log_file[ne_col].values
        
        # Calculate quantiles
        for q in np.arange(0.05, 1.05, 0.05):
            upper = np.quantile(rate, 1 - q/2)
            lower = np.quantile(rate, q/2)
            
            if q == 1:
                lower = lower + 0.03
                upper = upper - 0.03
            
            reassortment_data.append({
                'time': mrsi - pd.Timedelta(days=shift * 365),
                'quantile': q,
                'upper': upper,
                'lower': lower,
                'isconstant': is_constant,
                'name': 'both',
                'method': 'constant rho' if is_constant else 'dependent rho(t)',
                'alpha': 1.0 if q == 1 else 0.2
            })
    
    return pd.DataFrame(reassortment_data)


def plot_reassortment_rate(ax, reassortment_df, smoothed_cases, mrsi, lineage_colors):
    """Plot reassortment rate over time with confidence bands"""
    # Plot reassortment rate bands
    for q in reassortment_df['quantile'].unique():
        q_data = reassortment_df[reassortment_df['quantile'] == q]
        ax.fill_between(q_data['time'], q_data['lower'], q_data['upper'],
                       alpha=q_data['alpha'].iloc[0], 
                       color=lineage_colors.get(q_data['name'].iloc[0], '#4DAF4A'))
    
    # Add case data if provided
    if smoothed_cases is not None:
        for case_type in smoothed_cases['type'].unique():
            type_data = smoothed_cases[smoothed_cases['type'] == case_type]
            color = '#E41A1C' if case_type == 'HPAI' else '#377EB8'
            ax.plot(type_data['date'], np.log(type_data['positivity']), 
                   color=color, linewidth=0.5, label=case_type)
    
    # Format axes
    ax.set_xlim(pd.to_datetime("2021-09-01"), mrsi)
    ax.set_ylim(-4, np.log(5))
    
    # Set y-axis ticks and labels for log scale
    y_ticks = [np.log(0.05), np.log(0.1), np.log(0.2), np.log(0.4), 
               np.log(0.8), np.log(1.6), np.log(3.2)]
    y_labels = ["0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "3.2"]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels)
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Reassortment rate")
    ax.legend(loc='upper left')
    
    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='y')


def load_clade_events(path, is_constant=True):
    """Load HPAI clade events data"""
    events_file = f"{path}/combined/HLHxNx.constant.clades.tsv"
    
    if os.path.exists(events_file):
        events = pd.read_csv(events_file, sep='\t')
        # Convert height to time
        mrsi = pd.to_datetime("2025-02-17")
        events['Time'] = mrsi - pd.to_timedelta(events['Height'] * 365.25, unit='D')
        # Filter events after 2020
        events = events[events['Time'] >= pd.to_datetime("2020-01-01")]
        return events
    return None


def calculate_event_distribution(hpai_events, segment_order):
    """Calculate distribution of reassortment events"""
    if hpai_events is None:
        return None, None
    
    lineage_types = ["HPAI", "LPAI"]
    event_types = ["HPAI", "HPAI+LPAI", "LPAI"]
    
    # Event distribution
    distr_data = []
    for sample in hpai_events['Sample'].unique():
        for lineage in lineage_types:
            lineage_events = hpai_events[(hpai_events['Sample'] == sample) & 
                                        (hpai_events['Lineage'] == lineage)]
            for event in event_types:
                if lineage not in event:
                    continue
                n_events = len(lineage_events[lineage_events['Event'] == event])
                distr_data.append({
                    'Sample': sample,
                    'Lineage': lineage,
                    'Event': event,
                    'n_events': n_events
                })
    
    # Co-reassortment analysis
    co_rea_data = []
    for sample in hpai_events['Sample'].unique():
        for lineage in lineage_types:
            lineage_events = hpai_events[(hpai_events['Sample'] == sample) & 
                                        (hpai_events['Lineage'] == lineage)]
            for event in event_types:
                if lineage not in event:
                    continue
                event_data = lineage_events[lineage_events['Event'] == event]
                
                for j, seg in enumerate(segment_order[1:], 1):  # Skip HA
                    # Count occurrences of segment index in Segments column
                    count = event_data['Segments'].apply(lambda x: str(j-1) in str(x)).sum()
                    co_rea_data.append({
                        'Sample': sample,
                        'Lineage': lineage,
                        'Event': event,
                        'n_events': count,
                        'segment': seg
                    })
    
    return pd.DataFrame(distr_data), pd.DataFrame(co_rea_data)


def plot_co_reassortment(ax, co_rea_df, lineage_colors):
    """Plot co-reassortment events by segment"""
    if co_rea_df is None:
        ax.text(0.5, 0.5, "No co-reassortment data available", 
                ha='center', va='center', transform=ax.transAxes)
        return
    
    # Calculate quantiles for HPAI lineage
    hpai_data = co_rea_df[co_rea_df['Lineage'] == 'HPAI']
    
    summary = []
    for segment in hpai_data['segment'].unique():
        for event in hpai_data['Event'].unique():
            seg_event_data = hpai_data[(hpai_data['segment'] == segment) & 
                                      (hpai_data['Event'] == event)]
            if len(seg_event_data) > 0:
                summary.append({
                    'segment': segment,
                    'Event': event,
                    'mean': seg_event_data['n_events'].mean(),
                    'lower': seg_event_data['n_events'].quantile(0.025),
                    'upper': seg_event_data['n_events'].quantile(0.975)
                })
    
    summary_df = pd.DataFrame(summary)
    
    # Rename events for clarity
    summary_df['evname'] = summary_df['Event'].replace({
        'HPAI': 'HPAI',
        'HPAI+LPAI': 'LPAI'
    })
    
    # Plot
    segments = summary_df['segment'].unique()
    x_pos = np.arange(len(segments))
    
    for i, event_name in enumerate(['HPAI', 'LPAI']):
        event_data = summary_df[summary_df['evname'] == event_name]
        offset = -0.15 if i == 0 else 0.15
        
        means = [event_data[event_data['segment'] == seg]['mean'].values[0] 
                if seg in event_data['segment'].values else 0 
                for seg in segments]
        lowers = [event_data[event_data['segment'] == seg]['lower'].values[0] 
                 if seg in event_data['segment'].values else 0 
                 for seg in segments]
        uppers = [event_data[event_data['segment'] == seg]['upper'].values[0] 
                 if seg in event_data['segment'].values else 0 
                 for seg in segments]
        
        color = lineage_colors.get(event_name, '#4DAF4A')
        ax.errorbar(x_pos + offset, means, yerr=[np.array(means) - np.array(lowers), 
                                                  np.array(uppers) - np.array(means)],
                   fmt='o', color=color, label=f'from {event_name}', capsize=3)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(segments)
    ax.set_xlabel("Segment")
    ax.set_ylabel("Events")
    ax.legend(title="Lineage origin of segment")
    
    # Style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.3, axis='y')


def create_multi_segment_plot(path, segment_order, is_constant=True):
    """Create plot with all segments using baltic_bacter"""
    fig, axes = plt.subplots(2, 4, figsize=(12, 8))
    axes = axes.flatten()
    
    lineage_colors = {
        "HPAI": "#E41A1C",
        "LPAI": "#377EB8",
        "unknown": "#4DAF4A"
    }
    
    # Get mrsi_dec
    _, _, _, _, _, mrsi_dec, _, _, _, _, _ = define_constants()
    
    for i, segment in enumerate(segment_order[1:]):  # Skip HA
        ax = axes[i]
        ll = load_and_process_tree(path, segment, is_constant)
        draw_segment_tree(ax, ll, mrsi_dec, lineage_colors, segment, linewidth=0.5)
        ax.set_title(segment)
    
    # Add legend in the last subplot
    axes[-1].axis('off')
    # Create legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=lineage_colors['HPAI'], label='HPAI'),
                      Patch(facecolor=lineage_colors['LPAI'], label='LPAI'),
                      Patch(facecolor=lineage_colors['unknown'], label='Unknown')]
    axes[-1].legend(handles=legend_elements, loc='center', 
                   title='Pathogenicity\nbased on HA', fontsize=12)
    
    plt.tight_layout()
    return fig


def create_main_figure(path, is_constant=True):
    """Create the main composite figure (Figure 2) using baltic_bacter"""
    # Load all necessary data
    clades, rate_shifts, mrsi, mrsi_hpai, mrsi_lpai, mrsi_dec, segment_order, \
    methods_colors, species_colors, clade_colors, lineage_colors = define_constants()
    
    log_file = load_log_file(path, is_constant)
    cases = load_cases_data(path)
    smoothed_cases = smooth_case_data(cases)
    hpai_events = load_clade_events(path, is_constant)
    distr_df, co_rea_df = calculate_event_distribution(hpai_events, segment_order)
    
    # Create figure with subplots
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 3, width_ratios=[2, 1, 1], height_ratios=[1, 1])
    
    # Panel A: Main HA tree with clade coloring
    ax_tree = fig.add_subplot(gs[:, 0])
    ll_ha = load_and_process_tree(path, "HA", is_constant)
    draw_ha_tree_with_clades(ax_tree, ll_ha, mrsi_dec, linewidth=0.5)
    ax_tree.set_title("A", loc='left', fontweight='bold')
    
    # Panel B: Co-reassortment
    ax_corea = fig.add_subplot(gs[0, 1:])
    plot_co_reassortment(ax_corea, co_rea_df, lineage_colors)
    ax_corea.set_title("B", loc='left', fontweight='bold')
    
    # Panel C: NP tree
    ax_tree_c = fig.add_subplot(gs[1, 1])
    ll_np = load_and_process_tree(path, "NP", is_constant)
    draw_segment_tree(ax_tree_c, ll_np, mrsi_dec, lineage_colors, "NP", linewidth=0.5)
    ax_tree_c.set_title("C", loc='left', fontweight='bold')
    
    # Panel D: PB2 tree
    ax_tree_d = fig.add_subplot(gs[1, 2])
    ll_pb2 = load_and_process_tree(path, "PB2", is_constant)
    draw_segment_tree(ax_tree_d, ll_pb2, mrsi_dec, lineage_colors, "PB2", linewidth=0.5)
    ax_tree_d.set_title("D", loc='left', fontweight='bold')
    
    plt.tight_layout()
    return fig


def main():
    """Main execution function"""
    setup_matplotlib()
    
    # Set path - adjust as needed
    path = "/Users/nmueller/Documents/github/CoInfection-Material/Applications/H5N1NorthAmerica"
    
    # Process both constant and dependent models
    is_constant = True
    print(f"\nProcessing {'constant' if is_constant else 'dependent'} model...")
    
    # Load data
    clades, rate_shifts, mrsi, mrsi_hpai, mrsi_lpai, mrsi_dec, segment_order, \
    methods_colors, species_colors, clade_colors, lineage_colors = define_constants()
    
    try:
        log_file = load_log_file(path, is_constant)
        cases = load_cases_data(path)
        smoothed_cases = smooth_case_data(cases)
        reassortment_df = calculate_reassortment_rates(log_file, rate_shifts, mrsi, is_constant)
        
        # Create reassortment rate plot
        
        # Create main composite figure with baltic_bacter trees
        fig2 = create_main_figure(path, is_constant)
        fig2.savefig(f"{output_dir}Figure2.pdf", bbox_inches='tight')
        plt.show()
        plt.close(fig2)

        
        # # Create multi-segment plot with baltic_bacter trees
        # fig3 = create_multi_segment_plot(path, segment_order, is_constant)
        # fig3.savefig(f"{output_dir}h5n1_all_segment_trees_{'constant' if is_constant else 'dependent'}.pdf",
        #             bbox_inches='tight')
        # plt.close(fig3)
        
        # print(f"Successfully created plots for {'constant' if is_constant else 'dependent'} model")
        
    except Exception as e:
        print(f"Error processing {'constant' if is_constant else 'dependent'} model: {e}")
        import traceback
        traceback.print_exc()
        
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()