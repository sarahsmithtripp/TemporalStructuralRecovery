"""
This script reads in output files from pre and post values to investigate potential reasons for differences.

"""
import rasterio
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    ##Read in with numpy to see where differences are 
    xr_test = rasterio.open('xr_test.tif').read(1) #just read for the predisturbance values
    np_test = rasterio.open('np_test.tif').read(1)
    diff = xr_test[0, :,:] - np_test
    mask = diff < -5
    masked = np.where(mask, obs, np.nan) ## mask to the value 
    # with rasterio.open("diff.tif", 'w', **save) as dst:
    #     dst.write(diff, indexes = 1)

    #get GYC vals for these pixels 
    dist_year_mask = np.where(mask, dist_year, np.nan)
    # with rasterio.open("gyc.tif", 'w', **save) as dst:
    #     dst.write(dist_year_mask, indexes = 1)

    # Extract the time series for these pixels
    time, height, width = masked.shape
    df = pd.DataFrame(masked.reshape([time,-1]).T).dropna()
    df['id'] = df.index + 1
    df.set_index('id')
    # Define the columns to pivot on
    pivot_columns = df.melt(id_vars = 'id', value_name = 'measure')
    sample_set = np.random.choice(pivot_columns['id'].unique(), 
                                10, replace =False)
    sample = pivot_columns[pivot_columns['id'].isin(sample_set)]

    in_dat = sample.groupby('id')
    # Create a new figure and axis for the plot
    fig, ax = plt.subplots()

    # Iterate through groups and plot a line for each group
    for name, group in in_dat:
        ax.plot(group['variable'], group['measure'], label=name)

    # Add a legend to show the ID labels
    #ax.legend(title='ID')

    # Set the labels for the x-axis and y-axis
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # Show the plot
    plt.show()
