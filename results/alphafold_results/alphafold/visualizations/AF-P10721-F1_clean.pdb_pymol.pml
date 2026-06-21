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

load "data/AF-P10721-F1_clean.pdb", protein
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

load "data/AF-P10721-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [4944,5403,5404,5500,4950,4947,4954,5497,4403,4405,4406,4530,5499,5501,4358,4359,4360,4922,4790,4902,4897,4354,4353,4356,4357,4352,4943,4919,4923,4896,4542,4952,5498,5494,5496,4545,5473,5478,5480,4679,4710,4711,4788,5483,4056,4677,4881,4712] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.416,0.278,0.702]
select surf_pocket2, protein and id [6295,5748,5742,5746,6311,6291,6293,6297,6300,6308,5284,6304,5733,5770,5737,5739,4096,4097,5287,5283,5285,5286] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.878]
select surf_pocket3, protein and id [4381,4383,5504,5505,4379,4617,4619,4621,4563,5495,4031,4034,4546,4547,4651,4656,4680,4647,4648,4649,4650,4653,4620,5487,5502,5488,5490,5503,4040,4042,4032,4388,4544] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.380]
select surf_pocket4, protein and id [5183,5185,5194,5222,5437,5439,5440,6426,6453,6455,6458,5441,5452,5433,5151,5149,5157,5162,5152,6454,6427,6428] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.620,0.361]
select surf_pocket5, protein and id [2850,2851,2852,2539,2541,2544,2847,2849,2613,2802,2797,2879,2821,2791,2588,2589] 
set surface_color,  pcol5, surf_pocket5 
   

deselect

orient
