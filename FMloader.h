#pragma once

#ifndef LOADER_H
#define LOADER_H
#endif // !LOADER_H


#include <map>
#include <string>
#include <vector>
#include <filesystem>

/* 
    OpenVSP parser for DCS or other C++ applications
    by Loki_v2

    Licensed under GNU Public License v3
    You may use and modify this code as long as you give credit to contributors.

*/


/*
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
 */

/*
    MAIN TODOs:
    - find a way to include control deflection as a parameter (maybe through clever file naming?)
    - exception handling as new ways to break this are discovered
    - maybe extend the airframe class to make the element class to make things cleaner
*/


#define verbose false
#define silent false

// which params do we care about loading
#define interesting_params {"CDi","CDo","CDtot","CDtrefftz","CFx","CFy","CFz","CL","CMx","CMy","CMz","Cmx","Cmy","Cmz","CS","Cms", "E","FC_AoA_","FC_Beta_","FC_Bref_","FC_Cref_","FC_Mach_","FC_Pitch_Rate","FC_ReCref_","FC_Rho_","FC_Roll__Rate","FC_Sref_","FC_Vinf_","FC_Xcg_","FC_Yaw___Rate","FC_Ycg_","FC_Zcg_","L / D"}



    class FMDataLoader              // container of all known data
    {

    public:
        class FMElementData         // container of all data on a single element
        {
        public:
            std::string name;
            std::vector<double> ELdeflects;  //  list of known dimension values for the element
            std::vector<double> ELmachs;
            std::vector<double> ELalphas;
            std::vector<double> ELbetas;
            
            // element constructors and setters
            FMElementData::FMElementData(std::string);                                                                       // basic constructor that just takes the element name
            FMElementData::FMElementData(std::string n, std::string p, double d, double m, double a, double b, double val);  // full constructor that initializes one param value
            void insert(std::string param, double deflect, double mach, double alpha, double beta, double value);            // add a new param value to an element

            // element accessors
            double lookup(std::string param, double deflect, double mach, double alpha, double beta);                        // get one param value
            std::vector<double> lookupParamVector(std::string param);                                                        // get the complete ordered list of values for a single param

        private:
            std::map<std::tuple<std::string, double, double, double, double>, double> data_frame;
        };

        std::vector<double> AFdeflects;     // list of known dimension values for the airframe
        std::vector<double> AFmachs;
        std::vector<double> AFalphas;
        std::vector<double> AFbetas;

        // airframe total polars constructors and setters
        FMDataLoader::FMDataLoader();                                                                                         // basic default constructor
        FMDataLoader::FMDataLoader(std::filesystem::path p);                                                                  // constructor which collects all properly located and structured files and sets up the class's internal data frame, airframe_polars
        void insertPolar(std::string param, double deflect, double mach, double alpha, double beta, double value);            // add a new polar value at given defletion/mach/alpha/beta

        // convenient accessors for polars and individual values 
        double  FMDataLoader::getFMParam(std::string element, std::string param, double deflect, double mach, double alpha, double beta); // get a single parameter value for an element designated by name (as opposed to object reference)
        double FMDataLoader::getPolar(std::string param, double deflect, double mach, double alpha, double beta);             // get individual polar value

        // entire class instance exporters and util functions
        FMElementData& getCompleteElement(std::string element);                                                               // get one entire element object
        std::vector<double>  FMDataLoader::getFMParamVector(std::string element, std::string param);                          // get aero coeffiecient as flat std::vector of choice for given element from loaded CSV data
        std::vector<std::string> FMDataLoader::ListElementNames();




    private:

        std::vector<std::filesystem::path> foundfilelist; // list of processed files
        std::vector<FMElementData> elements;              // list of known elements
        std::map<std::tuple<std::string, double, double, double, double>, double> airframe_polars;  // main store of polars

        void FMDataLoader::loadVSPcsv(std::filesystem::path);  // utility function to load flight data straight from OpenVSP history.csv during construction
        void FMDataLoader::loadcsv(std::filesystem::path);  // utility function to load flight data straight from OpenVSP history.csv during construction
        static std::vector<double> FMDataLoader::getNearest(std::vector<double>& dim, double key);
       // std::vector<std::pair<std::string, std::vector<double>>> FMDataLoader::loadcsv(std::filesystem::path);  // TODO utility function to bulk load tabular CSV data
    };



