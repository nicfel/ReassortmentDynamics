#!/usr/bin/env python3
"""
Script to plot raw case data and compare it to smoothed predictors used in GLM analyses.

This script loads the APHIS surveillance data and extracts the smoothed predictors
from the XML files to create comparison plots showing how the raw data is processed
into the predictors used in the GLM analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import os
import re
from scipy import stats
import seaborn as sns

def setup_plotting():
    """Set up matplotlib defaults - matching plot_glm.py styling"""
    typeface = 'helvetica'
    plt.rcParams['font.weight'] = 300
    plt.rcParams['axes.labelweight'] = 300
    plt.rcParams['font.family'] = typeface
    plt.rcParams['font.size'] = 12  # Consistent base font size
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['savefig.bbox'] = 'tight'

def load_raw_case_data(path):
    """
    Load the raw APHIS surveillance data
    
    Args:
        path (str): Path to the project directory
        
    Returns:
        pd.DataFrame: Raw case data with processed dates
    """
    csv_path = os.path.join(path, 'tables', 'APHIS_WildBirdAvianInfluenzaSurveillanceDashboard_with_flyways.csv')
    
    if not os.path.exists(csv_path):
        # Try without flyways suffix
        csv_path = os.path.join(path, 'tables', 'APHIS_WildBirdAvianInfluenzaSurveillanceDashboard.csv')
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Could not find APHIS data file at {csv_path}")
    
    print(f"Loading raw case data from: {csv_path}")
    cases = pd.read_csv(csv_path)
    
    # Convert dates to datetime
    cases['date'] = pd.to_datetime(cases['Date_Collected'], format='%Y-%m-%d')
    cases['decimal_date'] = cases['date'].apply(lambda x: x.year + (x.dayofyear - 1) / 365.25)
    
    return cases

def extract_smoothed_predictors(path):
    """
    Extract smoothed predictor values from XML files
    
    Args:
        path (str): Path to the project directory
        
    Returns:
        dict: Dictionary containing predictor arrays for each type
    """
    xml_file_path = os.path.join(path, 'xmls', 'HPAI_HLHxNx.glm.rep0.xml')
    
    if not os.path.exists(xml_file_path):
        raise FileNotFoundError(f"Could not find XML file at {xml_file_path}")
    
    print(f"Loading smoothed predictors from: {xml_file_path}")
    
    with open(xml_file_path, 'r') as xml_file:
        xml_content = xml_file.read()
    
    # Extract predictor values using regex
    predictors = {}
    
    # Define predictor patterns
    predictor_patterns = {
        'lpai': r'<stateNode id="lpai" spec="parameter.RealParameter" value="([^"]+)"/>',
        'hpai': r'<stateNode id="hpai" spec="parameter.RealParameter" value="([^"]+)"/>',
        'total': r'<stateNode id="total" spec="parameter.RealParameter" value="([^"]+)"/>',
        'overlap': r'<stateNode id="overlap" spec="parameter.RealParameter" value="([^"]+)"/>',
        'lpai_nosummer': r'<stateNode id="lpai_nosummer" spec="parameter.RealParameter" value="([^"]+)"/>',
        'h5_lpai': r'<stateNode id="h5_lpai" spec="parameter.RealParameter" value="([^"]+)"/>'
    }
    
    # Extract rate shifts for time points
    rate_shifts_match = re.findall(r'<stateNode id="rateShifts" spec="RealParameter" value="([^"]+)"/>', xml_content)
    if rate_shifts_match:
        rate_shifts = np.array([float(x) for x in rate_shifts_match[0].split()])
    else:
        # Fallback rate shifts if not found in XML
        rate_shifts = np.array([0, 0.105936073059361, 0.211872146118721, 0.317808219178082, 0.423744292237443, 
                               0.529680365296804, 0.635616438356164, 0.741552511415525, 0.847488584474886, 
                               0.953424657534247, 1.05936073059361, 1.16529680365297, 1.27123287671233, 
                               1.37716894977169, 1.48310502283105, 1.58904109589041, 1.69497716894977, 
                               1.80091324200913, 1.90684931506849, 2.01278538812785, 2.11872146118721, 
                               2.22465753424658, 2.33059360730594, 2.4365296803653, 2.54246575342466, 
                               2.64840182648402, 2.75433789954338, 2.86027397260274, 2.9662100456621, 
                               3.07214611872146, 3.17808219178082, 3.28401826484018, 3.38995433789954, 
                               3.4958904109589, 3.60182648401827, 3.70776255707763, 3.81369863013699, 
                               3.91963470319635, 4.02557077625571, 4.13150684931507])
    
    # Extract predictor values
    for pred_name, pattern in predictor_patterns.items():
        matches = re.findall(pattern, xml_content)
        if matches:
            values = np.array([float(x) for x in matches[0].split()])
            # Standardize values (log transform and standardize as done in GLM)
            min_val = np.min(values[values > 0])
            log_values = np.log(values + min_val)
            log_values = (log_values - np.mean(log_values)) / np.std(log_values)
            predictors[pred_name] = {
                'raw_values': values,
                'standardized_values': log_values,
                'time_points': rate_shifts
            }
    
    # Calculate time points in calendar years (assuming max date around 2025)
    mrsi = 2025.12877  # Most recent sampling date
    predictors['time_points_calendar'] = mrsi - rate_shifts
    
    return predictors

def calculate_raw_case_metrics(cases, time_window_days=46):
    """
    Calculate raw case metrics similar to the smoothing process in buildXmls.R
    
    Args:
        cases (pd.DataFrame): Raw case data
        time_window_days (int): Time window for smoothing (default 46 days = 2 * min(diff(rate_shifts)))
        
    Returns:
        pd.DataFrame: Raw case metrics by date
    """
    print("Calculating raw case metrics...")
    
    # Get date range
    min_date = cases['date'].min()
    max_date = cases['date'].max()
    
    raw_metrics = []
    
    # Calculate metrics for each day with sliding window
    for current_date in pd.date_range(min_date + timedelta(days=time_window_days//2), 
                                     max_date - timedelta(days=time_window_days//2), 
                                     freq='D'):
        
        # Define time window
        window_start = current_date - timedelta(days=time_window_days//2)
        window_end = current_date + timedelta(days=time_window_days//2)
        
        # Get cases in window
        window_cases = cases[(cases['date'] >= window_start) & (cases['date'] <= window_end)]
        
        if len(window_cases) == 0:
            continue
        
        # Calculate metrics
        total_AIV = sum(window_cases['Final_IAV'] == 'Detected')
        total_H5 = sum(window_cases['Final_H5'] == 'Detected')
        total_HPAI = sum((window_cases['Final_H5'] == 'Detected') & 
                        (window_cases['Final_Pathogenicity'] == 'High Path AI'))
        
        # Calculate LPAI (H5 detected but not high path)
        total_LPAI = total_H5 - total_HPAI
        
        # Calculate prevalences
        total_samples = len(window_cases)
        AIV_prevalence = total_AIV / total_samples if total_samples > 0 else 0
        HPAI_prevalence = total_HPAI / total_samples if total_samples > 0 else 0
        LPAI_prevalence = total_LPAI / total_samples if total_samples > 0 else 0
        
        # Calculate overlap (LPAI prevalence * HPAI prevalence)
        overlap = LPAI_prevalence * HPAI_prevalence
        
        raw_metrics.append({
            'date': current_date,
            'decimal_date': current_date.year + (current_date.dayofyear - 1) / 365.25,
            'total_AIV': total_AIV,
            'total_H5': total_H5,
            'total_HPAI': total_HPAI,
            'total_LPAI': total_LPAI,
            'total_samples': total_samples,
            'AIV_prevalence': AIV_prevalence,
            'HPAI_prevalence': HPAI_prevalence,
            'LPAI_prevalence': LPAI_prevalence,
            'overlap': overlap
        })
    
    return pd.DataFrame(raw_metrics)

def create_comparison_plots(raw_metrics, predictors, output_path):
    """
    Create comparison plots showing raw vs smoothed predictors
    
    Args:
        raw_metrics (pd.DataFrame): Raw case metrics
        predictors (dict): Smoothed predictor data
        output_path (str): Path to save plots
    """
    print("Creating comparison plots...")
    
    # Create figure with subplots
    fig, axes = plt.subplots(3, 2, figsize=(15, 12))
    fig.suptitle('Raw Case Data vs Smoothed Predictors Used in GLM Analysis', fontsize=16, fontweight='bold')
    
    # Define colors - matching plot_glm.py color scheme
    colors = {
        'raw': '#2c3e50',      # Dark blue-grey (main branch color from plot_glm.py)
        'smoothed': '#C85A3C', # Strong terracotta (reassortment rate color)
        'LPAI': '#2ca02c',     # Green for LPAI (matches plot_glm.py)
        'HPAI': '#d62728',     # Red for HPAI (matches plot_glm.py)
        'Total': '#4E79A7'     # Blue (matches plot_glm.py)
    }
    
    # Plot 1: LPAI cases (log-standardized counts)
    ax1 = axes[0, 0]
    # Log transform the counts, adding 1 to avoid log(0), then standardize
    log_lpai_counts = np.log(raw_metrics['total_LPAI'] + 1)
    log_lpai_counts_std = (log_lpai_counts - np.mean(log_lpai_counts)) / np.std(log_lpai_counts)
    ax1.plot(raw_metrics['date'], log_lpai_counts_std, 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw LPAI (log-standardized counts)')
    ax1.set_title('LPAI Cases (Log-standardized Counts)', fontweight='bold')
    ax1.set_ylabel('Standardized Log(Cases + 1)')
    ax1.legend()
    
    # Plot 2: HPAI cases (log-standardized counts)
    ax2 = axes[0, 1]
    # Log transform the counts, adding 1 to avoid log(0), then standardize
    log_hpai_counts = np.log(raw_metrics['total_HPAI'] + 1)
    log_hpai_counts_std = (log_hpai_counts - np.mean(log_hpai_counts)) / np.std(log_hpai_counts)
    ax2.plot(raw_metrics['date'], log_hpai_counts_std, 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw HPAI (log-standardized counts)')
    ax2.set_title('HPAI Cases (Log-standardized Counts)', fontweight='bold')
    ax2.set_ylabel('Standardized Log(Cases + 1)')
    ax2.legend()
    
    # Plot 3: LPAI log-standardized prevalence comparison
    ax3 = axes[1, 0]
    # Log transform prevalence, adding small value to avoid log(0), then standardize
    log_lpai_prevalence = np.log(raw_metrics['LPAI_prevalence'] + 1e-6)
    log_lpai_prevalence_std = (log_lpai_prevalence - np.mean(log_lpai_prevalence)) / np.std(log_lpai_prevalence)
    ax3.plot(raw_metrics['date'], log_lpai_prevalence_std, 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw LPAI (log-standardized prevalence)')
    
    if 'lpai' in predictors:
        # Convert time points back to dates for plotting
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax3.plot(time_dates, predictors['lpai']['raw_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='Smoothed LPAI Predictor')
    
    ax3.set_title('LPAI Log-standardized Prevalence: Raw vs Smoothed', fontweight='bold')
    ax3.set_ylabel('Standardized Log(Prevalence + 1e-6)')
    ax3.legend()
    
    # Plot 4: HPAI log-standardized prevalence comparison
    ax4 = axes[1, 1]
    # Log transform prevalence, adding small value to avoid log(0), then standardize
    log_hpai_prevalence = np.log(raw_metrics['HPAI_prevalence'] + 1e-6)
    log_hpai_prevalence_std = (log_hpai_prevalence - np.mean(log_hpai_prevalence)) / np.std(log_hpai_prevalence)
    ax4.plot(raw_metrics['date'], log_hpai_prevalence_std, 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw HPAI (log-standardized prevalence)')
    
    if 'hpai' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax4.plot(time_dates, predictors['hpai']['raw_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='Smoothed HPAI Predictor')
    
    ax4.set_title('HPAI Log-standardized Prevalence: Raw vs Smoothed', fontweight='bold')
    ax4.set_ylabel('Standardized Log(Prevalence + 1e-6)')
    ax4.legend()
    
    # Plot 5: Total AIV log-standardized prevalence comparison
    ax5 = axes[2, 0]
    # Log transform prevalence, adding small value to avoid log(0), then standardize
    log_total_prevalence = np.log(raw_metrics['AIV_prevalence'] + 1e-6)
    log_total_prevalence_std = (log_total_prevalence - np.mean(log_total_prevalence)) / np.std(log_total_prevalence)
    ax5.plot(raw_metrics['date'], log_total_prevalence_std, 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw Total AIV (log-standardized prevalence)')
    
    if 'total' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax5.plot(time_dates, predictors['total']['raw_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='Smoothed Total Predictor')
    
    ax5.set_title('Total AIV Log-standardized Prevalence: Raw vs Smoothed', fontweight='bold')
    ax5.set_ylabel('Standardized Log(Prevalence + 1e-6)')
    ax5.legend()
    
    # Plot 6: Overlap log-standardized comparison
    ax6 = axes[2, 1]
    # Log transform overlap, adding small value to avoid log(0), then standardize
    log_overlap = np.log(raw_metrics['overlap'] + 1e-8)
    log_overlap_std = (log_overlap - np.mean(log_overlap)) / np.std(log_overlap)
    ax6.plot(raw_metrics['date'], log_overlap_std, 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw Overlap (log-standardized)')
    
    if 'overlap' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax6.plot(time_dates, predictors['overlap']['raw_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='Smoothed Overlap Predictor')
    
    ax6.set_title('Overlap Log-standardized: Raw vs Smoothed', fontweight='bold')
    ax6.set_ylabel('Standardized Log(Overlap + 1e-8)')
    ax6.legend()
    
    # Format x-axes and add timeline shading like plot_glm.py
    timewidth = 0.5  # Half-year intervals for shading
    fromval = 2021.5  # Start of data
    toval = 2025.0    # End of data range
    
    for ax in axes.flat:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
        ax.tick_params(axis='x', rotation=45)
        
        # Add timeline shading like plot_glm.py
        for i in np.arange(fromval, toval, 2 * timewidth):
            # Convert to datetime for shading
            start_date = pd.to_datetime(f'{int(i)}-01-01') + pd.Timedelta(days=(i % 1) * 365)
            end_date = start_date + pd.Timedelta(days=timewidth * 365)
            ax.axvspan(start_date, end_date, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)
        
        # Style spines like plot_glm.py
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.3)
        ax.set_axisbelow(True)
    
    # Set consistent y-axis ranges for comparison
    # Get the overall range from all log-transformed and standardized data
    all_log_values = np.concatenate([
        log_lpai_counts_std, log_hpai_counts_std, 
        log_lpai_prevalence_std, log_hpai_prevalence_std, log_total_prevalence_std, log_overlap_std
    ])
    
    # Also include predictor values if available
    if 'lpai' in predictors:
        all_log_values = np.concatenate([all_log_values, predictors['lpai']['raw_values']])
    if 'hpai' in predictors:
        all_log_values = np.concatenate([all_log_values, predictors['hpai']['raw_values']])
    if 'total' in predictors:
        all_log_values = np.concatenate([all_log_values, predictors['total']['raw_values']])
    if 'overlap' in predictors:
        all_log_values = np.concatenate([all_log_values, predictors['overlap']['raw_values']])
    
    # Remove any infinite or NaN values
    all_log_values = all_log_values[np.isfinite(all_log_values)]
    
    # Set consistent y-axis range with some padding
    y_min = np.min(all_log_values) - 0.5
    y_max = np.max(all_log_values) + 0.5
    
    # Apply consistent y-axis range to all subplots
    for ax in axes.flat:
        ax.set_ylim(y_min, y_max)
    
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_path, 'raw_vs_smoothed_predictors.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot to: {output_file}")
    
    plt.show()

def create_standardized_comparison_plots(raw_metrics, predictors, output_path):
    """
    Create plots showing standardized (log-transformed) predictors
    
    Args:
        raw_metrics (pd.DataFrame): Raw case metrics
        predictors (dict): Smoothed predictor data
        output_path (str): Path to save plots
    """
    print("Creating standardized comparison plots...")
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Standardized Predictors Used in GLM Analysis\n(Log-transformed and Standardized)', fontsize=14, fontweight='bold')
    
    # Standardize raw data for comparison
    raw_standardized = raw_metrics.copy()
    
    # LPAI standardization
    min_lpai = np.min(raw_metrics['LPAI_prevalence'][raw_metrics['LPAI_prevalence'] > 0])
    if min_lpai > 0:
        log_lpai = np.log(raw_metrics['LPAI_prevalence'] + min_lpai)
        raw_standardized['LPAI_std'] = (log_lpai - np.mean(log_lpai)) / np.std(log_lpai)
    
    # HPAI standardization
    min_hpai = np.min(raw_metrics['HPAI_prevalence'][raw_metrics['HPAI_prevalence'] > 0])
    if min_hpai > 0:
        log_hpai = np.log(raw_metrics['HPAI_prevalence'] + min_hpai)
        raw_standardized['HPAI_std'] = (log_hpai - np.mean(log_hpai)) / np.std(log_hpai)
    
    # Total standardization
    min_total = np.min(raw_metrics['AIV_prevalence'][raw_metrics['AIV_prevalence'] > 0])
    if min_total > 0:
        log_total = np.log(raw_metrics['AIV_prevalence'] + min_total)
        raw_standardized['Total_std'] = (log_total - np.mean(log_total)) / np.std(log_total)
    
    # Overlap standardization
    min_overlap = np.min(raw_metrics['overlap'][raw_metrics['overlap'] > 0])
    if min_overlap > 0:
        log_overlap = np.log(raw_metrics['overlap'] + min_overlap)
        raw_standardized['Overlap_std'] = (log_overlap - np.mean(log_overlap)) / np.std(log_overlap)
    
    # Plot standardized predictors - matching plot_glm.py color scheme
    colors = {'raw': '#2c3e50', 'smoothed': '#C85A3C'}
    
    # LPAI standardized
    ax1 = axes[0, 0]
    ax1.plot(raw_standardized['date'], raw_standardized['LPAI_std'], 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw LPAI (Standardized)')
    if 'lpai' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax1.plot(time_dates, predictors['lpai']['standardized_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='GLM LPAI Predictor')
    ax1.set_title('Standardized LPAI Predictor', fontweight='bold')
    ax1.set_ylabel('Standardized Value')
    ax1.legend()
    
    # HPAI standardized
    ax2 = axes[0, 1]
    ax2.plot(raw_standardized['date'], raw_standardized['HPAI_std'], 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw HPAI (Standardized)')
    if 'hpai' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax2.plot(time_dates, predictors['hpai']['standardized_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='GLM HPAI Predictor')
    ax2.set_title('Standardized HPAI Predictor', fontweight='bold')
    ax2.set_ylabel('Standardized Value')
    ax2.legend()
    
    # Total standardized
    ax3 = axes[1, 0]
    ax3.plot(raw_standardized['date'], raw_standardized['Total_std'], 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw Total AIV (Standardized)')
    if 'total' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax3.plot(time_dates, predictors['total']['standardized_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='GLM Total Predictor')
    ax3.set_title('Standardized Total AIV Predictor', fontweight='bold')
    ax3.set_ylabel('Standardized Value')
    ax3.legend()
    
    # Overlap standardized
    ax4 = axes[1, 1]
    ax4.plot(raw_standardized['date'], raw_standardized['Overlap_std'], 
             color=colors['raw'], linewidth=2, alpha=0.7, label='Raw Overlap (Standardized)')
    if 'overlap' in predictors:
        time_dates = pd.to_datetime((predictors['time_points_calendar'] - 1970) * 365.25, unit='D')
        ax4.plot(time_dates, predictors['overlap']['standardized_values'], 
                 color=colors['smoothed'], linewidth=2, alpha=0.8, label='GLM Overlap Predictor')
    ax4.set_title('Standardized Overlap Predictor', fontweight='bold')
    ax4.set_ylabel('Standardized Value')
    ax4.legend()
    
    # Format x-axes and add timeline shading like plot_glm.py
    timewidth = 0.5  # Half-year intervals for shading
    fromval = 2021.5  # Start of data
    toval = 2025.0    # End of data range
    
    for ax in axes.flat:
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
        ax.tick_params(axis='x', rotation=45)
        
        # Add timeline shading like plot_glm.py
        for i in np.arange(fromval, toval, 2 * timewidth):
            # Convert to datetime for shading
            start_date = pd.to_datetime(f'{int(i)}-01-01') + pd.Timedelta(days=(i % 1) * 365)
            end_date = start_date + pd.Timedelta(days=timewidth * 365)
            ax.axvspan(start_date, end_date, facecolor='#E8ECF0', edgecolor='none', alpha=0.85, zorder=0)
        
        # Style spines like plot_glm.py
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(0.5)
        ax.spines['bottom'].set_linewidth(0.5)
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.3)
        ax.set_axisbelow(True)
    
    # Set consistent y-axis ranges for standardized comparison
    # Get the overall range from all standardized data
    all_std_values = np.concatenate([
        raw_standardized['LPAI_std'].dropna(), 
        raw_standardized['HPAI_std'].dropna(),
        raw_standardized['Total_std'].dropna(), 
        raw_standardized['Overlap_std'].dropna()
    ])
    
    # Also include predictor standardized values if available
    if 'lpai' in predictors:
        all_std_values = np.concatenate([all_std_values, predictors['lpai']['standardized_values']])
    if 'hpai' in predictors:
        all_std_values = np.concatenate([all_std_values, predictors['hpai']['standardized_values']])
    if 'total' in predictors:
        all_std_values = np.concatenate([all_std_values, predictors['total']['standardized_values']])
    if 'overlap' in predictors:
        all_std_values = np.concatenate([all_std_values, predictors['overlap']['standardized_values']])
    
    # Remove any infinite or NaN values
    all_std_values = all_std_values[np.isfinite(all_std_values)]
    
    # Set consistent y-axis range with some padding
    std_y_min = np.min(all_std_values) - 0.5
    std_y_max = np.max(all_std_values) + 0.5
    
    # Apply consistent y-axis range to all subplots
    for ax in axes.flat:
        ax.set_ylim(std_y_min, std_y_max)
    
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_path, 'standardized_predictors_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved standardized comparison plot to: {output_file}")
    
    plt.show()

def main():
    """Main function to run the analysis"""
    # Set up paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = script_dir
    
    # Set up plotting
    setup_plotting()
    
    try:
        # Load raw case data
        raw_cases = load_raw_case_data(script_dir)
        print(f"Loaded {len(raw_cases)} raw case records")
        
        # Extract smoothed predictors
        predictors = extract_smoothed_predictors(script_dir)
        print(f"Extracted {len(predictors)} predictor types")
        
        # Calculate raw case metrics
        raw_metrics = calculate_raw_case_metrics(raw_cases)
        print(f"Calculated metrics for {len(raw_metrics)} time points")
        
        # Create comparison plots
        create_comparison_plots(raw_metrics, predictors, output_path)
        create_standardized_comparison_plots(raw_metrics, predictors, output_path)
        
        # Print summary statistics
        print("\n=== Summary Statistics ===")
        print(f"Date range: {raw_metrics['date'].min()} to {raw_metrics['date'].max()}")
        print(f"Total samples: {raw_metrics['total_samples'].sum():,}")
        print(f"Total LPAI cases: {raw_metrics['total_LPAI'].sum():,}")
        print(f"Total HPAI cases: {raw_metrics['total_HPAI'].sum():,}")
        print(f"Mean LPAI prevalence: {raw_metrics['LPAI_prevalence'].mean():.4f}")
        print(f"Mean HPAI prevalence: {raw_metrics['HPAI_prevalence'].mean():.4f}")
        
    except Exception as e:
        print(f"Error: {e}")
        raise

if __name__ == "__main__":
    main()
