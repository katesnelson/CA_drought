import numpy as np
import gdal 
from gdal import *
import osr

#http://www.nass.usda.gov/Research_and_Science/Cropland/sarsfaqs2.php#Section1_11.0

d = gdal.Open('/data/emily/WF/kate/cv_lulc.tif', GA_ReadOnly).ReadAsArray()
d = d.astype('int16')

#na
d[d == 255] = -9999.0

#LU change flag
m = np.ma.masked_where(d==-9999.0, d)

delta_lc = []

for i in range(d.shape[0]-1):
    delta = d[i+1,:,:] - d[i,:,:]
    delta[delta > 0] = 1
    delta[delta < 0] = 1
    delta_lc.append(delta)

lu_change = np.asarray(delta_lc)  #lu in present changed relative to previous time step
#note, 2015 not in here b/c no 2016 to compare to
   


#CONSTRUCTION OF AGGREGATE LULC CLASS

#low number
d[d==4] = 4  #sorghum (cat 4)
d[d==2] = 4  #cotton (cat 4)
d[d == 1] = 3 #corn  (cat 3)
d[d==6]=4 #sunflower (cat 4)


#1: barren and fallow
d[d == 131] = 1 #barren
d[d == 61] = 1  #fallow


#2: non-edible grasses
d[d==59] = 2  #sod
d[d == 58] = 2 #clover
d[d==37] = 2 #hay
d[d==36]= 2  #alfalfa

#3: grains
d[d==225] = 3 #winterwheat/corn
d[d==226] = 3 #oats/corn
d[d==230] = 3 #lettuce/wheat
d[d == 205] = 3  #tritical
d[d==3] = 3  #rice 
d[d== 24] = 3  #winter wheat
d[d==12]=3 #sweet corn
d[d==13] = 3 #popcorn
d[d==21] = 3 #barley
d[d==22] = 3 #wheat
d[d==23] = 3 #wheat
d[d==24] = 3 #wheat
d[d==27] = 3 #rye
d[d==28] = 3 #oats
d[d==234] = 3 #wheat/sorghum
d[d==235] = 3 #barley/sorghum
d[d==236] = 3 #wheat/surghum
d[d==237] = 3 #barley/corn
d[d==238] = 3 #winter wheat/cotton
d[d==29] = 3 #millet


#4: row crops and veggies
d[d==42]=4 #dry beans
d[d==43] = 4 #potatoes
d[d==33] = 4 #safflower
d[d==54] = 4 #tomatos
d[d==31] = 4 #canola
d[d==45] = 4 #sugarcane
d[d==14] = 4 #mint
d[d==57] = 4 #herbs
d[d==208] = 4 #garlic
d[d==46] = 4 #sweet potatos
d[d==49] = 4 #onions
d[d==50] = 4 #cucs
d[d==53] = 4 #peas
d[d==206] = 4 #carrots
d[d==207] = 4 #asparagus
d[d==214] = 4 #broc
d[d==216] = 4 #peppers
d[d==219] = 4 #greens
d[d==222]=4 #squash
d[d==227] = 4 #lettuce
d[d==229]= 4 #pumpkin
d[d==243] = 4 #cabbage
d[d==244] = 4 #cauliflower
d[d==245] = 4 #celery
d[d==246] = 4 #radish
d[d==248] = 4 #eggplant
d[d==47] = 4 #misc veg and fruit
d[d==41] = 4 #sugar beets
d[d==224] = 4 #vetch
d[d==231] = 4 #lettuce, cantaloupe
d[d==232] = 4 #lettuce/cotton
d[d==38] = 4 #camelina, brassica family
d[d==247] = 4 #turnips


#5:  fruits and nuts (trees)
d[d==69] = 5 #grapes
d[d==75] = 5  #almonds
d[d==76] = 5  #walnuts
d[d==211] = 5  #olives
d[d==66]=5  #cherries
d[d==67]=5  #peaches
d[d==68] = 5 #apples
d[d==71] = 5 #other tree
d[d==72] = 5 #citrus
d[d==74] = 5 #pecans
d[d==77] = 5 #pears
d[d==204] = 5 #pistachios
d[d==212]=5 #oranges
d[d==220]=5 #plums
d[d==223] = 5 #apricots
d[d==218] = 5 #nectarines
d[d==48] = 5 #watermelon
d[d==55] = 5 #caneberries (?)
d[d==209] = 5 #cantaloupe
d[d==213] = 5 #honeydew
d[d==217] = 5 #pomegranites
d[d==221] = 5 #strawberries
d[d==242] = 5 #blueberries
d[d==44] = 5 #other crops
d[d==210] = 5 #prunes

#6: uncultivated cover
d[d==143]=6  #mixed forest
d[d==152]=6 #shrub
d[d==190]=6 #woodland wetland
d[d==195]=6 #herb wetlands
d[d==176]=6 #grassland
d[d==92] = 6 #aquaculture
d[d==70] = 6 #christmas trees
d[d==87] = 6 #wetlands
d[d==63] = 6 #woodland

d[d > 6] = 6  #includes 111 (open water), 112 (perennial ice/snow), 121-124 (developed), 141-142 (forest)


#crop diversification index
huc_idx = np.unique(huc12)
lulc_idx = np.unique(lulc)

cd_im = np.zeros(shape = lulc.shape)

for i in range(lulc.shape[0]):
    lc = lulc[i,:,:]    

    for h in range(len(np.unique(huc12))-1):
        data = lc[huc12 == huc_idx[h]]
        cdt = 0
        for l in range(len(np.unique(lulc))):
            total_area = data.shape[0]
            lulc_count = len(data[data == lulc_idx[l]])
            if lulc_count == 0:
                out == 0
            else:
                out = -(lulc_count/total_area)*np.log(lulc_count/total_area)
            cdt += out
        cd_im[i,:,:][huc12 == huc_idx[h]] = np.repeat(cdt, len(data))
        
cd_im[cd_im == 0] = -9999.0



#percent land use category by HUC
huc12 = gdal.Open('/data/emily/WF/kate/cv_huc.tif').ReadAsArray()
huc_idx = np.unique(huc12)
 
y=8
data = d[y,:,:]

out = []

for i in range(1,7): 
      
    category = np.zeros(data.shape)
    
    for h in range(len(np.unique(huc12))-1):
        
        huc_data = data[huc12 == huc_idx[h]]
        total = huc_data.shape[0]
        subset = huc_data[huc_data == i]
        ss = subset.shape[0]
        perc = np.true_divide(ss, total)
        category[huc12 == huc_idx[h]] = perc
    
    out.append(category)

outa = np.asarray(out)

#test = np.sum(outa, axis=0)  #more or less all one


#write out tiffs
infile = '/data/emily/WF/kate/cv.tif'
im = gdal.Open(infile).ReadAsArray()

#tvp_im
for i in range(outa.shape[0]):
    arr = outa[i,:,:]
    nodatav = -9999.0
    data = gdal.Open(infile)

    [cols, rows] = arr.shape
    trans = data.GetGeoTransform()
    proj = data.GetProjection()
    outfile = '/data/emily/WF/kate/cv_perc_lulc/perc_y15_cat' + str(i+1) + '.tif'

    outdriver = gdal.GetDriverByName('GTiff')
    outdata = outdriver.Create(str(outfile), rows, cols, 1, gdal.GDT_Float32)

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    outdata.SetProjection(srs.ExportToWkt())
    outdata.SetGeoTransform(trans)
    outdata.SetProjection(proj)
     
    #write array to file
    outdata.GetRasterBand(1).WriteArray(arr)
    outdata.GetRasterBand(1).SetNoDataValue(nodatav)
    outdata = None

    test = gdal.Open(outfile, GA_ReadOnly).ReadAsArray()
