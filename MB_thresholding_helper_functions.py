import pandas as pd
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt

def compute_angles(x, y):
    # compute angles between the intensity of the mouse and rabbit antibodies
    angles = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] == 0 and y[i] == 0:
            print("both x and y are zero for i = " + str(i))
            angle = -9
        elif x[i] == 0 or y[i] == 0:
            if x[i] == 0:
                angle = 90
            if y[i] == 0:
                angle = 0
        else:
            angle = math.degrees(math.atan(y[i] / x[i]))

        angles[i] = angle
    return angles

def compute_slopes_and_intercepts(x_t, y_t, theta1, theta2):
    # for use when shifting origin in intensity plots
    m1 = math.tan(math.radians(theta1))
    m2 = math.tan(math.radians(theta2))
    c1 = y_t - m1 * x_t
    c2 = y_t - m2 * x_t
    return m1, m2, c1, c2

def plotter_initial(x, y, angles, intensity, bins):
    # initial plotting for decisions on cut-off points (thresholds) to decide blob category
    # categories = mouse, rabbit, both
    plt.figure(figsize=(20, 7), facecolor='w')
    xy_max = max(np.percentile(x, 99.9), np.percentile(y, 99.9))
    plt.subplot(121)
    plt.hist(angles, bins=bins)
    plt.title(intensity + ' angles')
    plt.xlabel('Angle')
    plt.ylabel('No. blobs')
    
    plt.subplot(122)
    plt.hist2d(x, y, bins=200, norm=mpl.colors.LogNorm(), range= [[0, xy_max], [0, xy_max]], cmap='plasma')
    cb = plt.colorbar()
    cb.set_label('Log(N)')
    plt.title(intensity + ' per blob')
    plt.xlabel(intensity +': FilteredTexRed')
    plt.ylabel(intensity + ': FilteredAtto647')

    plt.show()
    plt.close()
    
def plotter_initial_new_origin(x, y, x_t, y_t, intensity):
    # for shifting origin in intensity plots
    plt.figure(figsize=(9, 7), facecolor='w')
    xy_max = max(np.percentile(x, 99.9), np.percentile(y, 99.9))
    
    plt.hist2d(x, y, bins=200, norm=mpl.colors.LogNorm(), range= [[0, xy_max], [0, xy_max]], cmap='plasma')
    cb = plt.colorbar()
    cb.set_label('Log(N)')
    plt.title(intensity + ' per blob (x_t = ' + str(x_t) + ', y_t = ' + str(y_t) + ')')
    plt.xlabel(intensity +': FilteredTexRed')
    plt.ylabel(intensity + ': FilteredAtto647')
    plt.axhline(y=y_t, color='k', linestyle='dashed')
    plt.axvline(x=x_t, color='k', linestyle='dashed')

    plt.show()
    plt.close()

def plotter_thresholds(x, y, angles, intensity, bins, low, high, name, dir_path, expt_name):
    # plotting with cut-off points (thresholds) displayed
    plt.figure(figsize=(20, 7), facecolor='w')
    xy_max = max(np.percentile(x, 99.9), np.percentile(y, 99.9))
    plt.subplot(121)
    plt.hist(angles, bins=bins)
    plt.title(intensity + ' angles (low = ' + str(low) + ', high = ' + str(high) + ')')
    plt.xlabel('Angle')
    plt.ylabel('No. blobs')
    plt.axvline(low, color='k', linestyle='dashed', linewidth=1.5)
    plt.axvline(high, color='k', linestyle='dashed', linewidth=1.5)
    
    low_x = [0]
    low_y = [0]
    low_x.append(max(x))
    low_y.append(math.tan(math.radians(low)) * max(x))
    
    high_x = [0]
    high_y = [0]
    high_x.append(max(y) / math.tan(math.radians(high)))
    high_y.append(max(y))
    
    
    plt.subplot(122)
    plt.hist2d(x, y, bins=200, norm=mpl.colors.LogNorm(), range= [[0, xy_max], [0, xy_max]], cmap='plasma')
    cb = plt.colorbar()
    cb.set_label('Log(N)')
    plt.title(name + ' stained: ' + intensity + ' per blob')
    plt.xlabel(intensity +': FilteredTexRed')
    plt.ylabel(intensity + ': FilteredAtto647')
    plt.plot(low_x, low_y, color='k', linestyle='dashed', linewidth=1.5)
    plt.plot(high_x, high_y, color='k', linestyle='dashed', linewidth=1.5)
    plt.savefig(dir_path + expt_name + '_' + name + '_thresholded.pdf', format='pdf')

    plt.show()
    plt.close()
    
def plotter_thresholds_new_origin(x, y, intensity, x_t, y_t, m1, m2, c1, c2, theta1, theta2, name, dir_path, expt_name):
    # plotting with cut-off points (thresholds) displayed, for cases with new origin
    plt.figure(figsize=(9, 7), facecolor='w')
    xy_max = max(np.percentile(x, 99.9), np.percentile(y, 99.9))
    
    low_x = [0]
    low_y = [c1]
    low_x.append(xy_max)
    low_y.append(m1 * xy_max + c1)
    
    high_x = [0]
    high_y = [c2]
    high_x.append(xy_max)
    high_y.append(m2 * xy_max + c2)
    
    plt.hist2d(x, y, bins=200, norm=mpl.colors.LogNorm(), range= [[0, xy_max], [0, xy_max]], cmap='plasma')
    cb = plt.colorbar()
    cb.set_label('Log(N)')
    plt.title(name + ' stained: ' + intensity + ' per blob (theta1 = ' + str(theta1) + ', theta2 = ' + str(theta2) + ')')
    plt.xlabel(intensity +': FilteredTexRed')
    plt.ylabel(intensity + ': FilteredAtto647')
    plt.plot(low_x, low_y, color='k', linestyle='dashed', linewidth=1.5)
    plt.plot(high_x, high_y, color='k', linestyle='dashed', linewidth=1.5)
    # plt.axhline(y=y_t, color='gray', linestyle='dashed')
    # plt.axvline(x=x_t, color='gray', linestyle='dashed')
    plt.savefig(dir_path + expt_name + '_' + name + '_' + intensity + '_thresholded.pdf', format='pdf')

    plt.show()
    plt.close()

def df_maker_cells(angles, cells, condition, low, high):
    # making data frame for cells based on the cut-off points (thresholds) selected 
    identification = ['empty' for i in range(len(angles))]
    for i in range(len(identification)):
        if angles[i] <= low:
            identification[i] = 'TexRed'
        elif angles[i] >= high:
            identification[i] = 'Atto647'
        else:
            identification[i] = 'Both'
    
    data = {'ImageNumber': cells['ImageNumber'],
           'ParentCells': cells['Parent_Cells'],
           'BlobCategory': identification}
    
    blobs_in_cells_df = pd.DataFrame(data, columns=['ImageNumber', 'ParentCells', 'BlobCategory'])
    images = list(set(blobs_in_cells_df['ImageNumber']))
    im_no = []
    cell_no = []
    blob_count = []

    for i in range(len(images)):
        temp = blobs_in_cells_df[blobs_in_cells_df['ImageNumber']==images[i]]
        cells = list(set(temp['ParentCells']))
        for j in range(len(cells)):
            temp2 = temp[temp['ParentCells']==cells[j]]
            for k in range(3):
                im_no.append(images[i])
                cell_no.append(cells[j])
            blob_count.append(sum(temp2['BlobCategory']=='TexRed'))
            blob_count.append(sum(temp2['BlobCategory']=='Atto647'))
            blob_count.append(sum(temp2['BlobCategory']=='Both'))

    category = ['TexRed', 'Atto647', 'Both'] * int(len(im_no)/3)
    cond = [condition] * int(len(im_no))
    
    data = {'Condition': cond,
            'ImageNumber': im_no,
            'ParentCell': cell_no,
            'BlobCount': blob_count,
            'Category': category}
    
    df = pd.DataFrame(data, columns = ['Condition','ImageNumber', 'ParentCell', 'BlobCount', 'Category'])
    return df

def df_maker_cells_new_origin(x, y, x_t, y_t, m1, m2, c1, c2, cells, condition):
    identification = ['empty' for i in range(len(x))]
    for i in range(len(identification)):
        if (m1 * x[i] + c1 - y[i] > 0) and (x[i] - x_t > 0):
            identification[i] = 'TexRed'
        elif (y[i] - m2 * x[i] - c2 > 0) and (y[i] - y_t > 0):
            identification[i] = 'Atto647'
        elif (y[i] - m2 * x[i] - c2 <= 0) and (m1 * x[i] + c1 - y[i] <= 0) and (x[i] - x_t > 0) and (y[i] - y_t > 0):
            identification[i] = 'Both'
        else:
            identification[i] = 'Neither'

    to_keep = [i for i in range(len(identification)) if identification[i] != 'Neither']
    
    data = {'ImageNumber': cells['ImageNumber'],
           'ParentCells': cells['Parent_Cells'],
           'BlobCategory': identification}
    
    blobs_in_cells_df = pd.DataFrame(data, columns=['ImageNumber', 'ParentCells', 'BlobCategory'])
    
    #blobs_in_cells_df = blobs_in_cells_df.iloc[to_keep]
    
    images = list(set(blobs_in_cells_df['ImageNumber']))
    im_no = []
    cell_no = []
    blob_count = []

    for i in range(len(images)):
        temp = blobs_in_cells_df[blobs_in_cells_df['ImageNumber']==images[i]]
        cells = list(set(temp['ParentCells']))
        for j in range(len(cells)):
            temp2 = temp[temp['ParentCells']==cells[j]]
            for k in range(3):
                im_no.append(images[i])
                cell_no.append(cells[j])
            blob_count.append(sum(temp2['BlobCategory']=='TexRed'))
            blob_count.append(sum(temp2['BlobCategory']=='Atto647'))
            blob_count.append(sum(temp2['BlobCategory']=='Both'))

    category = ['TexRed', 'Atto647', 'Both'] * int(len(im_no)/3)
    cond = [condition] * int(len(im_no))
    
    data = {'Condition': cond,
            'ImageNumber': im_no,
            'ParentCell': cell_no,
            'BlobCount': blob_count,
            'Category': category}
    
    df = pd.DataFrame(data, columns = ['Condition','ImageNumber', 'ParentCell', 'BlobCount', 'Category'])
    return df

def df_maker_tissues(angles, cells, condition, low, high):
    # making data frame for tissues based on the cut-off points (thresholds) selected 
    identification = ['empty' for i in range(len(angles))]
    for i in range(len(identification)):
        if angles[i] <= low:
            identification[i] = 'TexRed'
        elif angles[i] >= high:
            identification[i] = 'Atto647'
        else:
            identification[i] = 'Both'
    
    data = {'ImageNumber': cells['ImageNumber'],
           'BlobCategory': identification}
    
    blobs_in_cells_df = pd.DataFrame(data, columns=['ImageNumber', 'BlobCategory'])
    images = list(set(blobs_in_cells_df['ImageNumber']))
    im_no = []
    blob_count = []

    for i in range(len(images)):
        for k in range(3):
            im_no.append(images[i])
        
        temp = blobs_in_cells_df[blobs_in_cells_df['ImageNumber']==images[i]]
        blob_count.append(sum(temp['BlobCategory']=='TexRed'))
        blob_count.append(sum(temp['BlobCategory']=='Atto647'))
        blob_count.append(sum(temp['BlobCategory']=='Both'))

    category = ['TexRed', 'Atto647', 'Both'] * int(len(im_no)/3)
    cond = [condition] * int(len(im_no))
    
    data = {'Condition': cond,
            'ImageNumber': im_no,
            'BlobCount': blob_count,
            'Category': category}
    
    df = pd.DataFrame(data, columns = ['Condition','ImageNumber', 'BlobCount', 'Category'])
    return df

def df_maker_tissues_new_origin(x, y, x_t, y_t, m1, m2, c1, c2, cells, condition):
    identification = ['empty' for i in range(len(x))]
    for i in range(len(identification)):
        if (m1 * x[i] + c1 - y[i] > 0) and (x[i] - x_t > 0):
            identification[i] = 'TexRed'
        elif (y[i] - m2 * x[i] - c2 > 0) and (y[i] - y_t > 0):
            identification[i] = 'Atto647'
        elif (y[i] - m2 * x[i] - c2 <= 0) and (m1 * x[i] + c1 - y[i] <= 0) and (x[i] - x_t > 0) and (y[i] - y_t > 0):
            identification[i] = 'Both'
        else:
            identification[i] = 'Neither'

    to_keep = [i for i in range(len(identification)) if identification[i] != 'Neither']
    
    data = {'ImageNumber': cells['ImageNumber'],
           'ParentCells': cells['Parent_Cells'],
           'BlobCategory': identification}
    
    blobs_in_cells_df = pd.DataFrame(data, columns=['ImageNumber', 'ParentCells', 'BlobCategory'])
    
    #blobs_in_cells_df = blobs_in_cells_df.iloc[to_keep]
    
    images = list(set(blobs_in_cells_df['ImageNumber']))
    im_no = []
    cell_no = []
    blob_count = []

    for i in range(len(images)):
        temp = blobs_in_cells_df[blobs_in_cells_df['ImageNumber']==images[i]]
        cells = list(set(temp['ParentCells']))
        for j in range(len(cells)):
            temp2 = temp[temp['ParentCells']==cells[j]]
            for k in range(3):
                im_no.append(images[i])
                cell_no.append(cells[j])
            blob_count.append(sum(temp2['BlobCategory']=='TexRed'))
            blob_count.append(sum(temp2['BlobCategory']=='Atto647'))
            blob_count.append(sum(temp2['BlobCategory']=='Both'))

    category = ['TexRed', 'Atto647', 'Both'] * int(len(im_no)/3)
    cond = [condition] * int(len(im_no))
    
    data = {'Condition': cond,
            'ImageNumber': im_no,
            'ParentCell': cell_no,
            'BlobCount': blob_count,
            'Category': category}
    
    df = pd.DataFrame(data, columns = ['Condition','ImageNumber', 'ParentCell', 'BlobCount', 'Category'])
    return df