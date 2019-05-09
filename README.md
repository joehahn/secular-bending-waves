## secular-waves

by Joe Hahn,<br />
jhahn@spacescience.org,<br />
20 April 2019<br />
git branch=master


### bending waves

Descend into folder bending_waves to find the IDL code for simulating the 
spiral bending waves that are described in Hahn, 2007, ApJ, 665, 856;
the code discussed there was developed in 2003 I think.
I myself no longer have IDL running at my fingertips...hopefully this 
still runs 16 years later on your more modern IDL install, fingers crossed.

To execute the simulation shown in Fig 4-6 of the 2007 paper, enter the outgoing_waves folder,
start IDL, and then

    %run ring_master.pro

to execute the simulation. Or so I think...its been a while. Then

    %run plot_system.pro

to view the simulation output, and

    %run paper_figs.pro

to generate the plots seen in Figs 3-6. I think.

And to execute the simulation seen in Figs 7-8, descend into the pan_in_out folder
and repeat the above. Good luck!


### density waves

Enter the density_waves folder for the IDL code that simulates the spiral density waves
described in Hahn, 2008, ApJ, 680, 1569, enter the outgoing_waves folder, and then

    %run ring_master.pro
    %run plot_system.pro
    %run paper_figs.pro

to generate the plots seen in Figs 3-6 of that paper. I think.


### kuiper belt

The above code was originally developed to study the waves that a young Neptune
may have launched in a more massive primordial Kuiper Belt. I haven't yet found that
simulation's initial conditions, but I did recover a nice animation of those
results, see files in the kuiper_belt folder


-jh
