from pymol import cmd,stored

set depth_cue, 1
set fog_start, 0.4

set_color b_col, [36,36,85]
set_color t_col, [10,10,10]
set bg_rgb_bottom, b_col
set bg_rgb_top, t_col      
set bg_gradient

set  spec_power  =  200
set  spec_refl   =  0

load "data/AF-Q92630-F1_clean.pdb", protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
#show_as cartoon, protein
show surface, protein
#set transparency, 0.15

show sticks, ligands
set stick_color, magenta




# SAS points

load "data/AF-Q92630-F1_clean.pdb_points.pdb.gz", points
hide nonbonded, points
show nb_spheres, points
set sphere_scale, 0.2, points
cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)


stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
lastSTP=stored.list[-1] # get the index of the last residue
hide lines, resn STP

cmd.select("rest", "resn STP and resi 0")

for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))



set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [1670,1682,1794,1798,1786,1789,1697,1698,1790,1299,636,1271,1274,1269,1275,1110,1112,1258,1248,1276,1669,1796,1797,949,1791,1249,673,693,641,642,644,671,648,692,800,635,639,816,817,818,1250,1252,1253,1281,1283,1293] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.365,0.278,0.702]
select surf_pocket2, protein and id [662,663,888,890,666,832,818,669,671,881,921,922,1811,1812,1813,950,885,1639,1640,1938,1819,1929,1913,1922,1817,1797,1641,1658,1943,1798] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.792,0.361,0.902]
select surf_pocket3, protein and id [1890,1886,1889,936,938,939,905,907,919,915,914,1897,1899,1901,1633,1818,1829,1831,1927,1824,1883,1867,940,1871,1817] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.533]
select surf_pocket4, protein and id [1975,1977,1962,1964,2220,2221,1316,1318,1320,1673] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.361]
select surf_pocket5, protein and id [230,233,250,1145,1147,753,1146,259,1159,1171,1162,575,555,556,1161,235,238,240,579,277] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.533,0.278]
select surf_pocket6, protein and id [2072,2082,1570,1572,1573,2076,2071,2073,2077,2078,2105,3038,2991,3010,3021,3026,3030] 
set surface_color,  pcol6, surf_pocket6 
   

deselect

orient
