#include "stdafx.h"
#include "FMLoader.h"
#include <stdexcept> 
#include <sstream>
#include <iostream>
#include <fstream>
//#include <utility> 
//#include <algorithm>

std::vector<std::string> paramlist interesting_params;

    // basic constructor for a single element
    FMDataLoader::FMElementData::FMElementData(std::string n) {
        name = n;
    }

    // full element constructor with a single param value
    FMDataLoader::FMElementData::FMElementData(std::string n, std::string p, double d, double m, double a, double b, double val) {
        name = n;
        insert(p, d, m, a, b, val);
    }

    // insert 1 new param value
    void FMDataLoader::FMElementData::insert(std::string param, double deflect, double mach, double alpha, double beta, double value) {

 
        if (find(deflects.begin(), deflects.end(), deflect) == deflects.end())  // if new element, add it to the unique list of elements
        {
            deflects.push_back(deflect);
        }

        if (find(machs.begin(), machs.end(), mach) == machs.end())  // if new mach, add it to the unique list of machs
        {
            machs.push_back(mach);
        }

        if (find(alphas.begin(), alphas.end(), alpha) == alphas.end())  // if new alpha, add it to the unique list of alphas
        {
            alphas.push_back(alpha);
        }

        if (find(betas.begin(), betas.end(), beta) == betas.end())  // if new beta, add it to the unique list of betas
        {
            betas.push_back(beta);
        }

        data_frame.insert({ make_tuple(param, deflect, mach, alpha, beta), value });
        if (verbose) printf("Inserted value %f for param %s in element %s\n", value, param.c_str(), name.c_str());
    }




    double FMDataLoader::FMElementData::lookup(std::string param, double deflect, double mach, double alpha, double beta) {
        // get a single parameter value at a given deflection, mach, alpha, beta
        // par names are same as OpenVSP nomenclature (CL, Cmx, CDi, etc)

        double result = 0.0;
        try {
            result = data_frame.at(make_tuple(param, deflect, mach, alpha, beta));
        }
        catch (const std::out_of_range& oor) {
            if (!silent) std::cout << "FM param not found for : " << name << ", " << param << ", " << deflect << ", " << mach << ", " << alpha << ", " << beta << std::endl;
        }
        return result;

    }




    std::vector<double> FMDataLoader::FMElementData::lookupParamVector(std::string param) {
        // get a well ordered std::vector of known values for the given parameter in this element
        // the default sort is betas then alphas then machs then deflections

        std::vector<double> result;
        std::tuple<std::string, double, double, double, double> key;

        for (auto& d : deflects)
        {
            for (auto& m : machs)
            {
                for (auto& a : alphas)
                {
                    for (auto& b : betas)
                    {
                        key = make_tuple(param, d, m, a, b);
                        try {
                            result.push_back(data_frame.at(key));
                        }
                        catch (const std::out_of_range& oor) {
                            if (!silent) std::cout << "FM param not found for : " << name << "," << param << "," << d << "," << m << "," << a << "," << b << std::endl;
                        }
                    }
                }
            }
        }
        return result;
    }

    FMDataLoader::FMDataLoader() {
        // scans the \Saved Games\DCS.openbeta\Mods\aircraft\Project-Lancaster\bin\ directory for csv files, assumed to be OpenVSP history files, and constructs an accessible data structure from them.

        std::string userdir = std::string(getenv("USERPROFILE"));
        std::filesystem::path p;

        if (!userdir.empty())
        {
            p = userdir + "\\Saved Games\\DCS.openbeta\\Mods\\aircraft\\" + MOD_FOLDER + "\\bin\\";
            //std::cout << "Loading FM data from ... " << p.string() << std::endl;
            if (verbose) printf("Loading FM data from %s\n", p.c_str());
        }
        else
        {
            if (!silent) printf("No user path");
        }

        for (const auto& entry : std::filesystem::directory_iterator(p))
        {

            if (entry.path().extension().compare(".csv") == 0)
            {
                foundfilelist.push_back(entry.path().stem());
                // printf("Loading %s \n", entry.path().string().c_str());

                FMDataLoader::loadVSPcsv(entry.path());
            }
        }
    }


    // mechanics of opening and reading a CSV into usable std::vectors
    // parts adapted from https ://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/ 
    std::vector<std::pair<std::string, std::vector<double>>> FMDataLoader::loadcsv(std::filesystem::path filename)
    {

        std::vector<std::pair<std::string, std::vector<double>>> result;

        std::string line, colname;
        double val;
        int num_lines = 0;

        std::ifstream myFile(filename.string());

        if (!myFile.is_open()) throw std::runtime_error("LancFM unable to load file parameter file");
        std::cout << "Processing.." << filename << std::endl;


        // Read the column names
        if (myFile.good())
        {
            //    std::cout << "file is good" << std::endl;
                // Extract the first line in the file
            std::getline(myFile, line);

            // Create a std::stringstream from line
            std::stringstream ss(line);

            // Extract each column name
            while (std::getline(ss, colname, ',')) {

                // Initialize and add <colname, int std::vector> pairs to result
                result.push_back({ colname,{} });
                //      std::cout << colname << std::endl;
            }
        }



        // Read data, line by line
        while (std::getline(myFile, line))
        {
            // Create a std::stringstream of the current line
            std::stringstream ss(line);

            // Keep track of the current column index
            int colIdx = 0;

            // Extract each integer
            while (ss >> val) {

                // Add the current integer to the 'colIdx' column's values std::vector
                result.at(colIdx).second.push_back(val);

                // If the next token is a comma, ignore it and move on
                if (ss.peek() == ',') ss.ignore();

                // Increment the column index
                colIdx++;
            }
            num_lines += 1;
        }

        if (!silent) std::cout << "Loaded " << num_lines << "lines" << std::endl;

        // Close file
        myFile.close();
        std::cout << "FM Loading Complete!";

        return result;
    }




    // load data directly from OpenVSP history.csv
    void FMDataLoader::loadVSPcsv(std::filesystem::path filename)
    {
        if(!silent) printf("Processing %s\n", filename.string().c_str());

        std::map<std::tuple< std::string, std::string, double, double, double, double>, double> result;
        // <parameter, map<key=<deflection, mach, alpha, beta>, value=parameter value>>>

        std::vector<std::pair<std::string, std::string>> pages; // break the file up into parseable page strings and keep track of which page is what type of VSP output
        int num_lines = 0;

        std::ifstream myFile(filename.string());

        if (!myFile.is_open()) throw std::runtime_error("LancFM unable to load file OpenVSP history file");

        // Read the column names
        if (myFile.good())
        {
            // Extract first column in the file so we can paginate it
            std::string line;
            int linecount = 0;

            std::string p = "";      // we'll store each page here temporarily
            std::string cur_page_name = "start of file";  // what type of VSP output is this page (elements? totals? other?)

            while (getline(myFile, line))       // first we'll break the text file into consumable Pages delimited by the VSP "Results_Name" line.
            {
                std::string r, t;               // temp string holders
                linecount++;

                std::stringstream ss(line); 
                std::getline(ss, r, ',');       // get first comma-separate token in line
              
                if (!r.compare("Results_Name")) // store the page start and corresponding VSP page name (e.g. VSPAERO_Comp_Load)
                {
                    if (p.length() > 0) pages.push_back(make_pair(p, cur_page_name));   // since we're now starting a new page, let's save the previous page to the vector of pages.  Note cur_page_name is still storing the previous pages name.
                    p = "";                                         // reset the temp page accumulator

                    std::getline(ss, cur_page_name, ',');           // grab the next csv token, which should be the value of Results_Name (e.g VSPAERO_Comp_Load)                            
                    if (verbose) printf("New VSP result page found starting line %d\n", linecount);
                }
                else // if the line isn't a result heading, append it to the page string for the last seen result
                {
                    p += line;
                    p += '\n';
                }
            }
            // the pages std::vector should now have pairs of strings holding the page contents and VSP name for each page.
            // the pages we care about are VSPAERO_Comp_Load, but let's keep references to all pages for future usecases
            // VSPAERO_Comp_Load pages have predictable structure that we can parse
            /*  Results_Name  -- std::stringstream starts here
                Results_Timestamp  -- dont care
                Results_Date       -- dont care
                Results_Time       -- dont care
                AoA
                Beta
                CDi
                CFx
                CFy
                CFz
                CL
                Cmx
                Cmy
                Cmz
                Comp_ID             -- dont care
                Comp_Name
                Cs
                Mach
             */

            std::vector<std::string> elem;      // remember the elements we've found in the results, this should be stable for any single VSP history file.

            for (const auto& page : pages)      // step through our std::vector of pages and process them into data
            {
                double deflection = 0.0; // TODO read this in somehow, its not obvious in the results
                double alpha;
                double beta;
                double mach;
                std::vector<double> deflects = { 0.0 };
                std::vector<double> alphas;
                std::vector<double> betas;
                std::vector<double> machs;
                std::string l,t;


                if (verbose) printf("Processing... %s\n", std::get<1>(page).c_str());

                // airframe-total parameter loading
                if (!std::get<1>(page).compare("VSPAERO_Polar")) {

                    if (verbose) printf("Found a useful page of polars %s \n", std::get<0>(page).c_str());
                    
                    std::vector<std::tuple<std::string, std::vector<double>>>  param_cache;

                    std::stringstream pagestream(std::get<0>(page));
                    {
                        std::stringstream linestream(l);
                        while (getline(linestream, t, ','))
                        {

                            if (t.compare("Alpha") == 0)
                            {
                                std::string x;
                                while (getline(linestream, x, ',')) // record alphas for this page in order 
                                {
                                    alphas.push_back(stod(x));
                                }
                            }
                            else if (!t.compare("Beta") == 0)
                            {
                                std::string x;
                                while (getline(linestream, x, ',')) // record betas for this page in order
                                {
                                    betas.push_back(stod(x));
                                }
                            }
                            else if (!t.compare("Mach") == 0)
                            {
                                std::string x;
                                while (getline(linestream, x, ',')) // record machs list for this page in order
                                {
                                    machs.push_back(stod(x));
                                }
                            }
                           // else if (t.compare("CDi") == 0 || t.compare("CDo") == 0 || t.compare("CDtot") == 0 || t.compare("CFx") == 0 || t.compare("CFy") == 0 || t.compare("CFz") == 0 || t.compare("CL") == 0 || t.compare("Cmx") == 0 || t.compare("Cmy") == 0 || t.compare("Cmz") == 0 || t.compare("Cms") == 0 || t.compare("CS"))
                            else if (std::find(paramlist.begin(), paramlist.end(), t) != paramlist.end())
                            {
                                std::string x;
                                std::vector<double> v;

                                while (getline(linestream, x, ','))
                                {
                                    param_cache.push_back(make_tuple((t), v));
                                }
                            }

                        }
                        // all lines vectorized, lets push them into final data frame
                        for (auto& param : param_cache)
                        {
                            int col_indx = 0;
                            for (auto& val : std::get<1>(param))
                            {
                                insertPolar(std::get<0>(param), deflection, machs.at(col_indx), alphas.at(col_indx), betas.at(col_indx), val);
                                col_indx++;
                            }
                        }
                    }
                }

                // element-wise parameter loading..
                if (!std::get<1>(page).compare("VSPAERO_Comp_Load"))
                {
                    if (verbose) printf("Found a useful result page %s \n", std::get<0>(page).c_str());

                    std::stringstream pagestream(std::get<0>(page));
                    while (getline(pagestream, l,'\n'))  // get the next line from the page
                    {
                        std::stringstream linestream(l);
                        while (getline(linestream, t, ','))
                        {
                            if (t.compare("AoA") == 0)
                            {
                                std::string x;
                                getline(linestream, x, ','); // record alpha for this page (its the same for all columns in page)
                                alpha = stod(x);
                            }
                            else if (t.compare("Beta") == 0)
                            {
                                std::string x;
                                getline(linestream, x, ','); // record beta for this page (its the same for all columns in page)
                                beta = stod(x);
                            }
                            else if (t.compare("Mach") == 0)
                            {
                                std::string x;
                                getline(linestream, x, ','); // record mach for this page (its the same for all columns in page)
                                mach = stod(x);
                            }
                            else if (t.compare("Comp_Name") == 0)
                            {
                                std::string x;
                                bool found = false;
                                while (getline(linestream, x, ','))
                                {
                                  //  printf("%s\t", x.c_str());
                                    found = false;
                                    for (auto& el : elem) {
                                        if (x.compare(el)==0)
                                        {
                                            found = true;
                                            //printf("Known elems %s",el.c_str());
                                        }
                                    }
                                    if (!found && x != "") {
                                        elem.push_back(x);  // hooray we identified a previously unknown element
                                        if (verbose) printf("New element identified %s\n", x.c_str());
                                    }
                                }
                            }
                            else if (std::find(paramlist.begin(), paramlist.end(), t) != paramlist.end())
                            //else if (t.compare("CDi") == 0 || t.compare("CFx") == 0 || t.compare("CFy") == 0 || t.compare("CFz") == 0 || t.compare("CL") == 0 || t.compare("Cmx") == 0 || t.compare("Cmy") == 0 || t.compare("Cmz") == 0 || t.compare("Cms") == 0)
                            {
                                std::string x;
                                int elem_iter = 0;
                                while (getline(linestream, x, ',')) 
                                {
                                    try
                                    {
                                        elements.push_back(FMElementData(elem.at(elem_iter), t, deflection, mach, alpha, beta, stod(x)));
                                        if (verbose) printf("Value %f pushed \n", stod(x));
                                        elem_iter++;
                                    }
                                    catch (const std::out_of_range& oor) {
                                        if (!silent) printf("Element not found for %d %f\n", elem_iter, stod(x));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        myFile.close();
        if (!silent) std::cout << "FM Loading Complete!";

    }

    FMDataLoader::FMElementData& FMDataLoader::getCompleteElement(std::string element)
    {
        int l = 0;
        int found_index = -1;

        for (auto& m : elements)
        {
            if (!(m.name.compare(element))) {
                found_index = l;
                break;
            }
            l++;
        }

        if (found_index < 0) { std::cout << element << " not found in data! Be sure to use same element names as the OpenVSP output"; }
        return elements.at(found_index);
    }


    std::vector<double> FMDataLoader::getFMParamVector(std::string element, std::string param) // get list of given parameter (Cl, Cd, etc) values for "element from CSV data
    {
        // We want to get back a Cl, Cd (CMy, etc) std::vector for each AeroElement (identified in the csv file)
        // The std::vector should be sorted first by deflection, then by mach, then by alpha, then (finally) by beta

        bool listuniques = false;
        std::vector<double> resvec;

        if (!(param.compare("AoA") || param.compare("Beta") || param.compare("Deflection") || param.compare("Mach")))
        {
            bool listuniques = true;
        }

        for (auto& m : elements)
        {
            if (!m.name.compare(element)) {
                resvec = m.lookupParamVector(param);
            }
        }

        return resvec;
    }


    double FMDataLoader::getFMParam(std::string element, std::string param, double deflect, double mach, double alpha, double beta) // get list of given parameter (Cl, Cd, etc) values for "element from CSV data
    {
        // We want to get back a Cl, Cd (CMy, etc) std::vector for each AeroElement (identified in the csv file)
        // The std::vector should be sorted first by deflection, then by mach, then by alpha, then (finally) by beta

        double result = 0.0;
        std::tuple<std::string, double, double, double, double> key = make_tuple(param, deflect, mach, alpha, beta);

        for (auto& m : elements)
        {
            if (!m.name.compare(element))
            {
                result = m.lookup(param, deflect, mach, alpha, beta);
            }
        }
        return result;
    }


    void FMDataLoader::insertPolar(std::string param, double deflect, double mach, double alpha, double beta, double value) 
    {
        try
        {
            airframe_polars.insert({ make_tuple(param, deflect, mach, alpha, beta), value });
            if (verbose) printf("Inserted value %f for param %s in airframe polars\n", value, param.c_str());
        }
        catch (const std::out_of_range& oor) {
            if (!silent) printf("Polar not found for %d %f\n", param, value);
        }
    }

    double FMDataLoader::getPolar(std::string param, double deflect, double mach, double alpha, double beta)
    {
        try {
            return airframe_polars.at(make_tuple(param, deflect, mach, alpha, beta));
        }
        catch (const std::out_of_range& oor) {
            if (!silent) printf("Polar value not found! %s %f %f %f %f\n", param.c_str(), deflect, mach, alpha, beta);
        }
        return 0.0;
    }
