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

load "data/AF-O14757-F1_clean.pdb", protein
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

load "data/AF-O14757-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1062,1063,1075,1180,1181,1182,715,716,288,289,290,418,1091,1090,694,706,690,673,693,709,273,667,668,710,707,708,1064,735,1176,1170,1172,1174,1175,1190,449,448,520,522,1178,1183,125,117,118,119,120,177,122,123,176,685,687,126,128,130,132,137,138,745,700,2201,2221,2218,2782,2779,2787,737,2765,2775,717,744,2236,2769,2774,2757,2767] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.416,0.278,0.702]
select surf_pocket2, protein and id [1050,1051,1062,1063,1075,1180,1181,1182,289,290,305,419,304,390,1178,1034,1198,1195,1196,1197,1174,132,138,149,177,153,155,360,383,1203,1204,1205] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.878]
select surf_pocket3, protein and id [269,547,549,551,257,653,655,657,659,546,548,656,658,660,44,45,68,69,70,71,72,73,75,213,576,580] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.380]
select surf_pocket4, protein and id [1362,1363,2799,1378,2599,2601,2603,1377,2565,2805,2809,2587,2602,2608,2798,1064,1065,1067,2801,1066,2800,2802,2803,1051,1049,1388] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.620,0.361]
select surf_pocket5, protein and id [2555,2869,2870,2871,2897,2899,2830,2566,2567,2569,2854,2853,2855,2862] 
set surface_color,  pcol5, surf_pocket5 
   

deselect

orient
