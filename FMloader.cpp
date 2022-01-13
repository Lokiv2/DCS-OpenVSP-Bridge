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

 
        if (find(ELdeflects.begin(), ELdeflects.end(), deflect) == ELdeflects.end())  // if new element, add it to the unique list of elements
        {
            ELdeflects.push_back(deflect);
        }

        if (find(ELmachs.begin(), ELmachs.end(), mach) == ELmachs.end())  // if new mach, add it to the unique list of machs
        {
            ELmachs.push_back(mach);
        }

        if (find(ELalphas.begin(), ELalphas.end(), alpha) == ELalphas.end())  // if new alpha, add it to the unique list of alphas
        {
            ELalphas.push_back(alpha);
        }

        if (find(ELbetas.begin(), ELbetas.end(), beta) == ELbetas.end())  // if new beta, add it to the unique list of betas
        {
            ELbetas.push_back(beta);
        }

        data_frame.insert({ make_tuple(param, deflect, mach, alpha, beta), value });
        if (verbose) printf("Inserted value %f for param %s in element %s\n", value, param.c_str(), name.c_str());
    }




    double FMDataLoader::FMElementData::lookup(std::string param, double deflect, double mach, double alpha, double beta) {
        // get a single parameter value at a given deflection, mach, alpha, beta
        // par names are same as OpenVSP nomenclature (CL, Cmx, CDi, etc)

        double result = 0.0;
        std::vector<double> x1;

        std::vector<double> def_proxy = getNearest(ELdeflects, alpha);
        std::vector<double> mach_proxy = getNearest(ELmachs, alpha);
        std::vector<double> alph_proxy = getNearest(ELalphas, alpha);
        std::vector<double> beta_proxy = getNearest(ELbetas, alpha);
        double a, b;

        try {
            a = data_frame.at(make_tuple(param, def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]));
            b = data_frame.at(make_tuple(param, def_proxy[1], mach_proxy[1], alph_proxy[1], beta_proxy[1]));
            result = (a + b) / 2.0;
        }
        catch (const std::out_of_range& oor) {
            if (!silent) printf("FM param not found for %s %s %f %f %f %f \n", name.c_str(), param.c_str(), deflect, mach, alpha, beta);
        }
        return result;

    }




    std::vector<double> FMDataLoader::FMElementData::lookupParamVector(std::string param) {
        // get a well ordered std::vector of known values for the given parameter in this element
        // the default sort is betas then alphas then machs then deflections

        std::vector<double> result;
        std::tuple<std::string, double, double, double, double> key;

        for (auto& d : ELdeflects)
        {
            for (auto& m : ELmachs)
            {
                for (auto& a : ELalphas)
                {
                    for (auto& b : ELbetas)
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

        AFdeflects.push_back(0.0); // temp until we have a way to get this from data
        std::vector<std::string> elem;      // remember the elements we've found in the results, this should be stable for any single VSP history file.
        std::vector<double> deflects = { 0.0 };
        std::vector<double> alphas;
        std::vector<double> betas;
        std::vector<double> machs;

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
            std::string cur_page_name = "VSPAERO_History";  // what type of VSP output is this page (elements? totals? other?)

            while (getline(myFile, line))       // first we'll break the text file into consumable Pages delimited by the VSP "Results_Name" line.
            {
                std::string r, t;               // temp string holders
                linecount++;

                std::stringstream ss(line); 
                std::getline(ss, r, ',');       // get first comma-separate token in line
              
                if (r.compare("Results_Name") == 0) // store the page start and corresponding VSP page name (e.g. VSPAERO_Comp_Load)
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



            for (const auto& page : pages)      // step through our std::vector of pages and process them into data
            {
                double deflection = 0.0; // TODO read this in somehow, its not obvious in the results
                double alpha;
                double beta;
                double mach;

                if (verbose) printf("Processing... \n%s\n", std::get<1>(page).c_str());

                // airframe-total parameter loading
                if (std::get<1>(page).compare("VSPAERO_Polar") == 0) 
                {
                    std::string l;
                    if (verbose) printf("Found a useful page of polars \n%s \n", std::get<0>(page).c_str());
                    std::vector<std::tuple<std::string, std::vector<double>>>  param_cache;

                    std::stringstream pagestream(std::get<0>(page));
                    while (getline(pagestream, l, '\n'))  // get the next line from the page
                    {
                        std::string t;
                        if (verbose) printf("\nline \t%s \n", l.c_str());

                        std::stringstream linestream(l);
                        getline(linestream, t, ',');

                        if (verbose) printf("Checking token \t%s", t.c_str());
                        if (t == "Alpha")
                        {
                            //   std::string x;
                            while (getline(linestream, t, ',')) // record alphas for this page in order 
                            {
                                if (std::find(AFalphas.begin(), AFalphas.end(), stod(t)) == AFalphas.end())
                                {
                                    AFalphas.push_back(stod(t));
                                    if (verbose)printf("adding new alphas to polars\t %f\n", stod(t));
                                }
                                alphas.push_back(stod(t));
                            }
                        }
                        else if (t == "Beta")
                        {
                            //    std::string x;
                            while (getline(linestream, t, ',')) // record betas for this page in order
                            {
                                if (std::find(AFbetas.begin(), AFbetas.end(), stod(t)) == AFbetas.end())
                                {
                                    AFbetas.push_back(stod(t));
                                    if (verbose)printf("adding new beta to polars\t %f\n", stod(t));

                                }
                                betas.push_back(stod(t));
                            }
                        }
                        else if (t == "Mach")
                        {
                            // std::string x;
                            while (getline(linestream, t, ',')) // record machs list for this page in order
                            {
                                if (std::find(AFmachs.begin(), AFmachs.end(), stod(t)) == AFmachs.end())
                                {
                                    AFmachs.push_back(stod(t));
                                    if(verbose) printf("adding new mach to polars\t %f\n", stod(t));

                                }
                                machs.push_back(stod(t));
                            }
                        }
                        // else if (t.compare("CDi") == 0 || t.compare("CDo") == 0 || t.compare("CDtot") == 0 || t.compare("CFx") == 0 || t.compare("CFy") == 0 || t.compare("CFz") == 0 || t.compare("CL") == 0 || t.compare("Cmx") == 0 || t.compare("Cmy") == 0 || t.compare("Cmz") == 0 || t.compare("Cms") == 0 || t.compare("CS"))
                        else if (std::find(paramlist.begin(), paramlist.end(), t) != paramlist.end())
                        {
                            std::string x;
                            std::vector<double> v;

                            while (getline(linestream, x, ','))
                            {
                                v.push_back(stod(x));
                            }
                            param_cache.push_back(make_tuple(t, v));
                        }
                    }

                    
                    // all lines vectorized, lets push them into final data frame
                    for (auto& param : param_cache)
                    {
                        int col_indx = 0;
                        for (auto& val : std::get<1>(param))
                        {
                        //    printf("trying to insert new polar %s %f", std::get<0>(param).c_str(), val);
                            insertPolar(std::get<0>(param), deflection, machs[col_indx], alphas[col_indx], betas[col_indx], val);
                            col_indx++;
                        }
                    }
                }

                // element-wise parameter loading..
                if (std::get<1>(page).compare("VSPAERO_Comp_Load")==0)
                {
                    std::string l;
                    if (verbose) printf("Found a useful result page %s \n", std::get<0>(page).c_str());

                    std::stringstream pagestream(std::get<0>(page));
                    while (getline(pagestream, l,'\n'))  // get the next line from the page
                    {
                        std::string t;
                        std::stringstream linestream(l);
                        while (getline(linestream, t, ','))
                        {
                            if (t.compare("AoA") == 0)
                            {
                                std::string x;
                                getline(linestream, x, ','); // record alpha for this page (its the same for all columns in page)
                                alpha = stod(x);
                                if (std::find(alphas.begin(), alphas.end(), alpha) != alphas.end())
                                {
                                    alphas.push_back(alpha);
                                }
                            }
                            else if (t.compare("Beta") == 0)
                            {
                                std::string x;
                                getline(linestream, x, ','); // record beta for this page (its the same for all columns in page)
                                beta = stod(x);
                                if (std::find(betas.begin(), betas.end(), beta) != betas.end())
                                {
                                    betas.push_back(beta);
                                }
                            }
                            else if (t.compare("Mach") == 0)
                            {
                                std::string x;
                                getline(linestream, x, ','); // record mach for this page (its the same for all columns in page)
                                mach = stod(x);
                                if (std::find(machs.begin(), machs.end(), mach) != machs.end())
                                {
                                     machs.push_back(mach);
                                }
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
        std::sort(AFdeflects.begin(), AFdeflects.end());
        std::sort(AFmachs.begin(), AFmachs.end());
        std::sort(AFalphas.begin(), AFalphas.end());
        std::sort(AFbetas.begin(), AFbetas.end());

        if (!silent) 
        {
            if (!silent) printf("FM Loading Complete!\n");
            if (!silent) printf("\nDeflections: \n");
            for (auto& x : AFdeflects) { printf("%f\t", x); }
            if (!silent) printf("\nMachs: \n");
            for (auto& x : AFmachs) { printf("%f\t", x); }
            if (!silent) printf("\nAlphas: \n");
            for (auto& x : AFalphas) { printf("%f\t", x); }
            if (!silent) printf("\nBetas: \n");
            for (auto& x : AFbetas) { printf("%f\t", x); }
            if (!silent) printf("\n");
        }

        for (auto& x : airframe_polars)
        {
            if (!silent) printf("%s, %f, %f, %f, %f = \t%f\n", std::get<0>(std::get<0>(x)), std::get<1>(std::get<0>(x)), std::get<2>(std::get<0>(x)), std::get<3>(std::get<0>(x)), std::get<4>(std::get<0>(x)), std::get<1>(x));
        }
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

        if (found_index < 0) { printf("\n not found in data! Be sure to use same element names as the OpenVSP output\n"); }
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
            if (m.name.compare(element)==0)
            {
                result = m.lookup(param, deflect, mach, alpha, beta);
            }
        }
        return result;
    }


    void FMDataLoader::insertPolar(std::string param, double deflect, double mach, double alpha, double beta, double value) 
    {
  //      printf("\nTrying to insert %s, %f, %f, %f, %f, %f\n", param, deflect, mach, alpha, beta, value);

        try
        {
            if (!std::get<1>(airframe_polars.insert({ make_tuple(param, deflect, mach, alpha, beta), value })))
            {
                airframe_polars[make_tuple(param, deflect, mach, alpha, beta)] = value;
                if (verbose) printf("value was already present, overwritten\n");
            }
            if(verbose) printf("\nInserted new value %f for param %s in airframe polars\n", value, param.c_str());
        }
        catch (const std::exception& e) { 
            if (!silent) printf("%s",e.what());
        }
    }

    double FMDataLoader::getPolar(std::string param, double deflect, double mach, double alpha, double beta)
    {
        std::vector<double> def_proxy = getNearest(AFdeflects, deflect);
        std::vector<double> mach_proxy = getNearest(AFmachs, mach);
        std::vector<double> alph_proxy = getNearest(AFalphas, alpha);
        std::vector<double> beta_proxy = getNearest(AFbetas, beta);

        double a, b;
        double result = 0.0;

        try {
            a = airframe_polars.at(make_tuple(param, def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]));
            if (verbose) printf("a test %s %f %f %f %f\n", param.c_str(), def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]);

            b = airframe_polars.at(make_tuple(param, def_proxy[1], mach_proxy[1], alph_proxy[1], beta_proxy[1]));
            if (verbose) printf("b test %s %f %f %f %f\n", param.c_str(), def_proxy[1], mach_proxy[1], alph_proxy[1], beta_proxy[1]);

            {
                if (verbose) printf("\nNearest are %f %f\n", a, b);
                result = (a + b) / 2.0;         // TODO proper interpolation
            }

        }
        catch (const std::out_of_range& oor) {
            if (!silent) printf("Polar value not found! %s %f %f %f %f\n", param.c_str(), deflect, mach, alpha, beta);
        }
        return result;
    }


    std::vector<double> FMDataLoader::getNearest(std::vector<double>& dim, double key)
    {
        std::vector<double>::iterator search_start = dim.begin();
        std::vector<double>::iterator search_end = dim.end();
        std::vector<double>::iterator midpoint;

        if (dim.size() < 1)
        {
            if (verbose) printf("getNearest: input array is empty!\n");
            return std::vector(0.0, 0.0);
        }

        if (verbose) printf("getting nearest %f in %d\n", key, dim.size() );

        std::vector<double> result{ *search_start, *search_start };  // initialize the return with the case where the dim has only 1 value

        // binary search through the list to find the nearest value match
        while (std::distance(search_start, search_end) > 1)
        {
            double min = *std::min_element(search_start, search_end);  // get the min and max values
            double max = *std::max_element(search_start, search_end);

            // catch boundary conditions
            if (key >= max) {
                result = { *search_end, *search_end };  // key is beyond max range of dim
        //        printf("off the max");
                break;
            }
            else if (key <= min) {
                result = { *search_start, *search_start }; //key is beyond min range of dim
          //      printf("off the min");
                break;
            }

            result =  {*search_start, *search_end};

            // binary search bifurcation
            midpoint = search_start + floor(std::distance(search_start, search_end) / 2.0);
            if (midpoint != search_start && midpoint != search_end)
            {
                if (key > *midpoint) search_start = midpoint;
                else search_end = midpoint;
            }
        } 

        if(verbose) printf("Got nearest %f %f\n", result[0], result[1]);
        return result;
    }
    