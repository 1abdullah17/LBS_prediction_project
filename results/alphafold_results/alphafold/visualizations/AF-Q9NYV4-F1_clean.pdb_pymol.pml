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

load "data/AF-Q9NYV4-F1_clean.pdb", protein
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

load "data/AF-Q9NYV4-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1592,1594,1597,1504,1494,1495,1496,1135,2932,2935,1599,1600,1601,1467,1483,1484,2940,443,434,489,490,436,439,2941,431,2944,449,452,458,459,466,468,448,474,1519,1520,1105,1102,1106,2930,2929,1077,1079,1080,605,602,604,606,607,865,1073,1075,1076,589,864] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.365,0.278,0.702]
select surf_pocket2, protein and id [3660,3685,3659,980,982,984,986,959,957,960,3686,3709,3710,3711,3674,3684,327,328,329,906,1038,983,985,987,1040,1041,1042,3702,763,729,733,734,3677,3678,641,3637,3638,688,689,690,3652,687,638,762] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.792,0.361,0.902]
select surf_pocket3, protein and id [2683,2685,2687,2688,2669,2678,2680,2682,2697,1915,1918,1921,1922,1925,1926,1927,1929,1394,1386,1388,1387,1399,2702,1361,1362,1419,1657,1924] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.533]
select surf_pocket4, protein and id [3390,3357,3359,3301,3426,3393,3395,3187,3272,3268,3280,3397,3274,3279] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.361]
select surf_pocket5, protein and id [2739,2601,2602,2603,2570,2737,2742,2705,2706,2707,2597,2576,2713] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.533,0.278]
select surf_pocket6, protein and id [2621,2625,2627,2632,2633,2253,2252,2212,2254,2626,2442,2443,2446] 
set surface_color,  pcol6, surf_pocket6 
   

deselect

orient
