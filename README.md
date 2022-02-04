# DSC-OpenVSP-Bridge
OpenVSP result file parser


    OpenVSP parser for DCS or other C++ applications
    by Loki_v2

    Licensed under GNU Public License v3
    You may use and modify this code as long as you give credit to contributors.



    CRASH COURSE:
    - compile FMloaded.h and FMloader.cpp into your project
    - place 1 or more unmodified OpenVSP history.csv files (can be renamed) in /config (or other path specificed in entry.lua) within your mod folder
    - AND/OR place 1 or more OpenVSP .polar files in /config.  When using .polar, which are random-number-of-spaces-separated, convert to proper comma-separated first.
    
    - In your FM code add:

    static FMDataLoader FMData()
    
    
    void ed_fm_configure(const char* cfg_path)
    {...
	FMdata = FMDataLoader(cfg_path); // load and store FM data
     ...}

    void ed_fm_simulate(double dt)
    {...

    	double Cd_actual = FMdata.getPolar("CL", 0.0, mach,  alpha_DEG, beta_DEG);
		double Cl_actual = FMdata.getPolar("CDtot", 0.0, mach, alpha_DEG, beta_DEG);
     ...}
     
    
     
     OR, if you want to have parameters for a given part of the airframe (AeroElement), you can do
     
     std::vector(FMDataLoader::AeroElement) aero;
     
     void ed_fm_configure(const char* cfg_path)
     {
     	...
		FMdata = FMDataLoader(cfg_path); // load and store FM data
		aero.push_back(AeroElement(FMdata, "Stabilizer", Vec3{ -12.5, 0.0, 2.6 }, Vec3{ 0.0, 0.0, 1.0 }, 3.7925, 0.9, -15.0, 20.0));  // configure the right stab, as an example
	...
      }
  
  and then 
  
 	void ed_fm_simulate(double dt)
	{
		for (AeroElement& elem : aero)
		{

		if (elem.name == "Stabilizer")
		{
			elem.deflect(elevator_DEG);
		}

		elem.aeroUpdate(Lancaster::mach, alpha_DEG, beta_DEG);   

		eDrag = elem.getDrag(Lancaster::ambientDensity_KgPerM3, Lancaster::totalVelocity_msec);
		eLift = elem.getLift(Lancaster::ambientDensity_KgPerM3, Lancaster::totalVelocity_msec);
		eMom = elem.getMoment(Lancaster::dynamicPressure_Nm2);
		eFcenter = elem.getForceCenter();

		add_force(eDrag, eFcenter);    // places the force at the appropriate airframe location. You'll need to figure something out to combine the force totals for DCS.
		add_force(eLift, eFcenter);    
		add_moment(eMom);

	}
     Enjoy your parameters!
    - if you want to be fancy you can grab the FMElementData objects and use their accessors within your aero element style code



    MAIN TODOs:
    - find a way to include control deflection as a parameter (maybe through clever file naming?)
    - exception handling as new ways to break this are discovered
    - maybe extend the airframe class to make the element class to make things cleaner
