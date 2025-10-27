import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from cycler import cycler

from matplotlib_custom_settings import *

xlsx = pd.ExcelFile(os.path.join(os.getcwd(), 'exported_results_new.xlsx'))

sheet_dict = {}
for sheet_name in xlsx.sheet_names:
    sheet_dict[sheet_name] = pd.read_excel(xlsx, sheet_name=sheet_name)

#%%

data_dict = {}

for sheet_name in xlsx.sheet_names:
    
    case_idx = re.search(r'case_(\d+)_', sheet_name).group(1)
    
    if case_idx[0] != '6':
    # if not any(x in case_idx[0] for x in ['1', '2', '6']):
    
        if '_2d' in sheet_name:
            
            # df = sheet_dict['case_4_2d']
            df = sheet_dict[sheet_name]
            
            # if not np.isnan(df["indep_var_1"].to_numpy()[0]):
                
            x_unique = np.sort(df["indep_var_1"].unique())
            y_unique = np.sort(df["indep_var_2"].unique())
            
            # Create 2D coordinate grids
            X, Y = np.meshgrid(x_unique, y_unique, indexing="ij")
            
            # Now reshape each variable into a 2D grid (matching the shape of X and Y)
            CDS_grid = df.pivot(index="indep_var_1", columns="indep_var_2", values="CDS").values
            CDS_ref_grid = df.pivot(index="indep_var_1", columns="indep_var_2", values="CDS_ref").values
            WMTO_grid = df.pivot(index="indep_var_1", columns="indep_var_2", values="WMTO").values
            WMTO_ref_grid = df.pivot(index="indep_var_1", columns="indep_var_2", values="WMTO_ref").values
            
            data_dict[case_idx] = {
                'CDS': CDS_grid,
                'CDS_ref': CDS_ref_grid,
                'WMTO': WMTO_grid,
                'WMTO_ref': WMTO_ref_grid,
            }
        
        elif '_1d' in sheet_name:
            
            df = sheet_dict[sheet_name]
            CDS_sequence = df['CDS'].to_numpy()
            CDS_ref_sequence = df['CDS_ref'].to_numpy()
            WMTO_sequence = df['WMTO'].to_numpy()
            WMTO_ref_sequence = df['WMTO_ref'].to_numpy()
            
            data_dict[case_idx] = {
                'CDS': CDS_sequence,
                'CDS_ref': CDS_ref_sequence,
                'WMTO': WMTO_sequence,
                'WMTO_ref': WMTO_ref_sequence,
            }

sys.exit()

#%%

plt.rcParams['axes.prop_cycle'] = cycler('color', ['#206095', '#a8bd3a', '#871a5b', '#f66068', '#05341A', '#27a0cc', '#003c57', '#22d0b6', '#746cb1', '#A09FA0'])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, ax = plt.subplots(figsize=(8, 6))

i = 0

for (key, value) in data_dict.items():
    
    color = colors[i]
    
    if value['CDS'].ndim == 2:

        x_unique[-1] -= 1e-6
        cs = ax.contour(
            value['WMTO'] / value['WMTO_ref'], value['CDS'] / value['CDS_ref'], X, levels=x_unique, colors=color, linewidths=0.8, zorder=10, extend='max', label=key
        )
        # ax.clabel(cs, fmt='%0.3f', fontsize=8)
        
        y_unique = y_unique.astype(float)
        y_unique[-1] -= 1e-6
        cs2 = ax.contour(
            value['WMTO'] / value['WMTO_ref'], value['CDS'] / value['CDS_ref'], Y, levels=y_unique, colors=color, linestyles='dashed', linewidths=0.8, zorder=10, extend='max',
        )
        # ax.clabel(cs2, fmt='%0.3f', fontsize=8)
        
        ax.scatter(value['WMTO'] / value['WMTO_ref'], value['CDS'] / value['CDS_ref'], color=color, marker='.', label=key)
        
        i += 1
    
    elif value['CDS'].ndim == 1:

        ax.plot(value['WMTO'] / value['WMTO_ref'], value['CDS'] / value['CDS_ref'], color=color, label=key)
        
        i += 1
    
        # ax.scatter(value['WMTO'] / value['WMTO_ref'], value['CDS'] / value['CDS_ref'], color=color, marker='.')

ax.set_xlabel("Net thrust (kN)")
ax.set_ylabel("Specific fuel consumption (g/kN/s)")
ax.spines[['right', 'top']].set_visible(False)
ax.tick_params(axis='y', which='both', right=False, length=0)
ax.tick_params(axis='x', which='both', length=0)

ax.set_xlim(left=1)
ax.set_ylim(bottom=1)

# proxy_gasturb = Patch(facecolor='blue', edgecolor='none', alpha=0.2, label='Flight envelope GasTurb')
# proxy_bada = Patch(facecolor='red',  edgecolor='none', alpha=0.2, label='Flight envelope BADA')
# handles, labels = ax.get_legend_handles_labels()
# handles.extend([proxy_gasturb, proxy_bada])
# ax.legend(handles=handles, labels=[*labels, 'Flight envelope GasTurb', 'Flight envelope BADA'], frameon=False)
ax.legend(frameon=False)

# ax.set_xscale('log')
# ax.set_yscale('log')

ax.set_aspect('equal')

plt.tight_layout()
plt.show()