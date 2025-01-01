
## Author: Junmin Wang
## Date: December 31st, 2024

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df_wide = pd.read_csv("/path/to/accuracies_comparison.csv")

# Reshape from wide to long format
df_long = df_wide.melt(var_name = 'Group', value_name = 'Value')
                       
# Calculate means and standard errors for error bars
means = df_long.groupby('Group')['Value'].mean()
std_errors = df_long.groupby('Group')['Value'].std()

# Create the plot
plt.figure(figsize=(6, 4))
sns.stripplot(x='Group', y='Value', data=df_long, jitter=True, color='grey', alpha=0.3, size=4)

# Add error bars (mean +/- std)
colors = ['blue', 'red']
for i, group in enumerate(df_long['Group'].unique()):
    plt.errorbar(x=i, y=means[group], yerr=std_errors[group], color=colors[i], capsize=10, label='_nolegend_')
    plt.errorbar(x=i, y=means[group], yerr=0, color=colors[i], capsize=10, label='_nolegend_')
    
# Add a bar
x1, x2 = 0, 1
y_max = max(df_long['Value']) + 0.1
plt.plot([x1, x2], [y_max, y_max], color='black', linewidth=1.5)

# Add significance annotation
p_value = 0.027
plt.text((x1 + x2) / 2, y_max + 0.02, f"p = {p_value:.3f}", ha='center', va='bottom', fontsize=10)

# Customize the plot
plt.xlabel('')
plt.ylabel('Accuracy')
plt.xticks(ticks=[0, 1], labels=['No Adjustment', 'ComBat-met'])
plt.xlim(-0.5, 1.5)
plt.ylim(0.4, 1.2)
plt.tight_layout()

# Save the plot
plt.savefig("/path/to/accuracies_jitter.png", dpi=300, bbox_inches='tight')
