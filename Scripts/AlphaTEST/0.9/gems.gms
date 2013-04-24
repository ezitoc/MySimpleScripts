(((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))
(                        INPUT SYSTEM FOR GeMS PROGRAM                         )
(                       comentaries are between brackets                       )
(                          caps can be lower or upper                          )
(                   caps must be upper for potential names                     )
(((((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))))

dimension 3
seed_first 130

(Agregamos moléculas)
>< read initial.xyz
sys add

(Grupos)
(Todo el sistema)
> sys  
  group 1 add
    
(Lo que se mueve)
> zrange 1 28.706020
  group 2 add
 
(Las dos capas fijas)
> group 1 
- group 2
  group 3 add

(Opciones)
> sys
set pbc   F F F

time step 0.020d0  (integration timestep [ps])

neigh verlet

interact 1 under tb

boost 1


(Velocidades y cuales se mueven)
> group 2
  set tempgdist 300  
  add ermak 300 1

(Selecciono la salida)
> group 1
  outfile 1 name traj.xyz
  outfile 1 pos 1
  outfile 1 each 50

  outfile 2 name ener.dat
  outfile 2 energy 1
  outfile 2 each 10

  outfile 3 name bias.dat
  outfile 3 bias
  outfile 3 each 10

  outfile 4 name temp.dat
  outfile 4 time temp 1
  outfile 4 each 10

> group 2
hyperd lp 3.75e8 15.9658671982 0.9
