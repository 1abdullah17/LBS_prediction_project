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

load "data/AF-Q9UKE5-F1_clean.pdb", protein
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

load "data/AF-Q9UKE5-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [165,160,163,164,214,216,217,159,162,774,775,780,161,1273,1280,1281,755,1270,316,776,760,1195,1272,791,792,795,1196,785,798,800,822,1169,1181,332,333,167,169,171,175,176] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.302,0.278,0.702]
select surf_pocket2, protein and id [4336,4792,4798,4790,4793,4815,4353,4324,4332,4334,4342,4031,4033,6243,4026,6237,6013,6235,6259,6245,6247,6249,4819,6261,6254,4040,4352,5574,5580,5201,5202,5572,5576,5565,5913,5915,5920,5907,4812,4813,4814,5579,5207,5932,5590,5928,5931,5924,5926,5558,5191,4327,4780] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.631,0.361,0.902]
select surf_pocket3, protein and id [4573,4576,4611,4072,4085,4069,4179,4181,4058,4068,4060,4061,4062,4063,4064,4066,4571,4603,4599,4572,4602,4164,4165] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.678,0.278,0.702]
select surf_pocket4, protein and id [788,791,792,790,818,854,829,785,2227,2229,2231,2232,2230,2259,2258,2293,2297,2256,2260,2261,2263,2222] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.682]
select surf_pocket5, protein and id [4265,4557,4514,4711,4556,4595,4665,4266,4294,4297] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.341]
select surf_pocket6, protein and id [6170,6173,6161,6163,3917,6229,3918] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.522,0.361]
select surf_pocket7, protein and id [597,598,599,297,299,600,601,603,259,261,285,286,296,113,114,625,260,263,102,262] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.596,0.278]
select surf_pocket8, protein and id [1274,1275,1282,754,1268,1284,1293,1294,1290,1289,453,484,486] 
set surface_color,  pcol8, surf_pocket8 
   

deselect

orient
