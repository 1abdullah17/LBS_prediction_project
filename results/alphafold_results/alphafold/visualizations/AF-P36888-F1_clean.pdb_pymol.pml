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

load "data/AF-P36888-F1_clean.pdb", protein
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

load "data/AF-P36888-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5049,5043,5046,5636,5635,5022,5013,5015,4445,5018,5607,5726,5728,5729,5757,4441,5722,5724,5727,4492,4494,4444,4619,5053,4442,4440,4446,4447,5745,5746,5754] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.278,0.702]
select surf_pocket2, protein and id [1712,1713,1711,1714,958,1511,1514,1518,1519,1520,994,1524,1539,1733,1685,1709,1719,1730,1734,2998,961,965,976,1720,992,970,974,985,959,1721,2606,2608,2604,1664,1666,1671,1665,2612,1542] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.576,0.361,0.902]
select surf_pocket3, protein and id [5860,5862,5864,4065,5867,5818,5821,5772,5868,4074,4062,4474,4475,5753,5769,5773,5774,4463,5779,5780,4113,5882,4096,4097,4114,4115,4116,4101,5584,4088,5731] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.616,0.278,0.702]
select surf_pocket4, protein and id [4871,4124,5699,5701,5706,5708,5711,4801,4870,4875,4876,4877,4802,4994,4996,5729,5715,5714,4768,5703,4123] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.792]
select surf_pocket5, protein and id [4221,4121,4126,4128,4130,4134,4099,4100,4098,4132,4106,4758,4230,4248,4249,4250,4253,4760,4232,4735,4737,4102,4105,4738,4793,4239] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.447]
select surf_pocket6, protein and id [1929,1934,1935,1936,1927,1937,1957,2189,2139,1889,1885] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.361]
select surf_pocket7, protein and id [449,718,617,618,619,632,443,446,492,488,476,477,624,690] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.447,0.278]
select surf_pocket8, protein and id [3088,3089,3260,3262,3574,3576,3578,3553,3522,3525,3554] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.792,0.361]
select surf_pocket9, protein and id [4190,4193,4164,5555,5542,5539,5538,5540,4145,6014,6020,6017,6021,5509,5512] 
set surface_color,  pcol9, surf_pocket9 
   

deselect

orient
