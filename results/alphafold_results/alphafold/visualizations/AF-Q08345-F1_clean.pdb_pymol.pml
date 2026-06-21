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

load "data/AF-Q08345-F1_clean.pdb", protein
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

load "data/AF-Q08345-F1_clean.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [3896,3899,3903,3926,3927,3928,4331,3178,3180,3867,3893,3871,3892,3179,3183,4431,4303,4465,4467,3846,4419,4420,4421,4422,4423,4425,4426,3464,4330,3185,3237,3239,3240,3182,3851,4317,4429] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.298,0.702]
select surf_pocket2, protein and id [3677,4240,3644,3651,4241,4242,3620,3614,3618,3653,3652,4408,4411,4259,4399,4405,3717,3725,4415,3645,4247,4254,4256,4251,4255] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.533,0.361,0.902]
select surf_pocket3, protein and id [227,228,229,230,231,232,234,235,393,1201,1203,403,404,1205,1206,405,386,388,1218,261,359,1210,682,628,631,683] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.565,0.278,0.702]
select surf_pocket4, protein and id [4414,4427,4428,4429,4430,3204,3198,3201,4413,3580,3584,3557,3217,3219,3213,4412,4422,3223,3197,3224,3186,3237,3240,3480,3481,3586,3621,3478,3588,3592,3594,3221,3222] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.878]
select surf_pocket5, protein and id [3488,3489,3808,3809,3823,3515,3249,3248,3125,3122,3123,3264,3141,3138,3096,3097,3095] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.533]
select surf_pocket6, protein and id [4740,4742,5232,5231,5233,5234,5236,4571,4579,4604,4770,4649,4947,4976,4978,4975,5000,5002,5003,4948,4887,4888,5091,5092,5235,5086,4973] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.490]
select surf_pocket7, protein and id [4844,4857,4859,4861,4862,4485,4459,4484,4486,3919,4464,4466,3918,4598,4807,4863,4865,4867,4868,4877,4869] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.329,0.278]
select surf_pocket8, protein and id [3852,3854,3853,3873,3850,4393,4394,3848,3855,3859,3697,3699,3737,3385,3386,3388,3700,4343,4346,4347,4348,4374,4377,4349,4351,3874,4340,4376] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.620,0.361]
select surf_pocket9, protein and id [212,216,221,268,271,219,196,194,199,203,266,211,189] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.631,0.278]
select surf_pocket10, protein and id [4138,5358,4113,4114,5357,5392,5393,4131,4133] 
set surface_color,  pcol10, surf_pocket10 
   

deselect

orient
