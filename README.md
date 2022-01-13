# DSC-OpenVSP-Bridge
OpenVSP result file parser


    OpenVSP parser for DCS or other C++ applications
    by Loki_v2

    Licensed under GNU Public License v3
    You may use and modify this code as long as you give credit to contributors.



    CRASH COURSE:
    - change the #def's below to match your folder and params needed
    - compile FMloaded.h and FMloader.cpp into your project
    - place 1 or more unmodified OpenVSP history.csv files (can be renamed) in /bin within your mod folder (usually alongside your EFM .dll)
    - In your FM code add:

    void fm_initialize()  // to initialize the data object and ingest OpenVSP files
    {...
        FMDataLoader FMdata = FMDataLoader(); // load and store FM data
    ...}

    void ed_fm_simulate(double dt)
    {...

    	double Cd_actual = FMdata.getPolar("CL", 0.0, mach,  alpha_DEG, beta_DEG);
		double Cl_actual = FMdata.getPolar("CDtot", 0.0, mach, alpha_DEG, beta_DEG);
     ...}

     Enjoy your parameters!
    - if you want to be fancy you can grab the FMElementData objects and use their accessors within your aero element style code



    MAIN TODOs:
    - find a way to include control deflection as a parameter (maybe through clever file naming?)
    - exception handling as new ways to break this are discovered
    - maybe extend the airframe class to make the element class to make things cleaner
