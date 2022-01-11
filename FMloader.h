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
    - you can now get individual FM params by calling 

        double FMDataLoader::getFMParam(std::string element, std::string param, double deflect, double mach, double alpha, double beta);            // get a single parameter value for an element designated by name (as opposed to object reference)
        double FMDataLoader::getPolar(std::string param, double deflect, double mach, double alpha, double beta);                                   // get individual airframe polar value

    - if you want to be fancy you can grab the FMElementData objects and use their accessors within your aero element style code
 */

/*
    MAIN TODOs:
    - find a way to include control deflection as a parameter (maybe through clever file naming?)
    - exception handling as new ways to break this are discovered
    - evaluate performance of the access operations and see if they need to be faster
*/



#define MOD_FOLDER "Project-Lancaster"  //TODO detect mod folder name within DCS.openbeta/mods/aircraft

#define verbose false
#define silent false

// which params do we care about loading
#define interesting_params std::vector<std::string>({"CDi","CDo","CDtot","CDtrefftz","CFx","CFy","CFz","CL","CMx","CMy","CMz","Cmx","Cmy","Cmz","CS","Cms", "E","FC_AoA_","FC_Beta_","FC_Bref_","FC_Cref_","FC_Mach_","FC_Pitch_Rate","FC_ReCref_","FC_Rho_","FC_Roll__Rate","FC_Sref_","FC_Vinf_","FC_Xcg_","FC_Yaw___Rate","FC_Ycg_","FC_Zcg_","L / D"})



    class FMDataLoader              // container of all known data
    {
    public:
        class FMElementData         // container of all known data on a single element
        {
        public:
            std::string name;
            std::vector<double> deflects;
            std::vector<double> machs;
            std::vector<double> alphas;
            std::vector<double> betas;
            
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

        // airframe total polars constructors and setters
        FMDataLoader::FMDataLoader();                                                                                         // constructor which collects all properly located and structured files and sets up the class's internal data frame, airframe_polars
        void insertPolar(std::string param, double deflect, double mach, double alpha, double beta, double value);            // add a new polar value at given defletion/mach/alpha/beta

        // convenient accessors for polars and individual values 
        double  FMDataLoader::getFMParam(std::string element, std::string param, double deflect, double mach, double alpha, double beta); // get a single parameter value for an element designated by name (as opposed to object reference)
        double FMDataLoader::getPolar(std::string param, double deflect, double mach, double alpha, double beta);             // get individual polar value
        std::vector<double> FMDataLoader::getPolarList(std::string param, double deflect, double mach, double alpha, double beta);     // get a given parameter as a flat vector



        // entire class instance exporter for use elsewhere
        FMElementData& getCompleteElement(std::string element);                                                               // get one entire element object
        std::vector<double>  FMDataLoader::getFMParamVector(std::string element, std::string param);                          // get aero coeffiecient as flat std::vector of choice for given element from loaded CSV data



    private:

        std::vector<std::filesystem::path> foundfilelist; // list of processed files
        std::vector<FMElementData> elements;              // list of known elements
        std::map<std::tuple<std::string, double, double, double, double>, double> airframe_polars;  // main store of polars

        void FMDataLoader::loadVSPcsv(std::filesystem::path);  // utility function to load flight data straight from OpenVSP history.csv during construction

        std::vector<std::pair<std::string, std::vector<double>>> FMDataLoader::loadcsv(std::filesystem::path);  // TODO utility function to bulk load tabular CSV data

    };

