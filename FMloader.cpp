#include "stdafx.h"
#include "FMLoader.h"
#include <stdexcept> 
#include <sstream>
#include <iostream>
#include <fstream>


std::vector<std::string> paramlist interesting_params;

// null constructor
    FMDataLoader::FMDataLoader() {
    }

// main container constructor
    FMDataLoader::FMDataLoader(std::filesystem::path p, bool elementwise) {
        // scans the /config/ directory for csv files, assumed to be OpenVSP history files, and constructs an accessible data structure from them.

        if (p == "")
        {
            if (!silent) printf("No user path");
        }

        if (!silent) printf("\nInitializing FM from %s\n", p.c_str());

        for (const auto& entry : std::filesystem::directory_iterator(p))
        {
            if (entry.path().extension().compare(".csv") == 0 && elementwise)
            {
                foundfilelist.push_back(entry.path().stem());
                FMDataLoader::loadVSPcsv(entry.path());
            }

            else if (entry.path().extension().compare(".polar") == 0 && !elementwise)
            {
                foundfilelist.push_back(entry.path().stem());
                FMDataLoader::loadcsv(entry.path());
            }
        }
    }


// basic empty constructor for a single element
    FMDataLoader::FMElementData::FMElementData(std::string n) {
        name = n;
    }

// full element constructor with a single param value
    FMDataLoader::FMElementData::FMElementData(std::string n, std::string p, double d, double m, double a, double b, double val) {
        name = n;
        insert(p, d, m, a, b, val);
    }



// Loaders

    void                FMDataLoader::loadcsv(std::filesystem::path filename)
    {
        // Loader of arbitrary csv (VSPAero .polar or hand-made data table
        // assumes VSP naming for parameters
        std::vector<std::pair<std::tuple<std::string, double, double, double>, double>> results;
        std::tuple<std::string, double, double, double> temp_key;

        std::string header, line;
        int num_lines = 0;
        int alpha_idx = 0;
        int beta_idx = 0;
        int mach_idx = 0;

        std::vector<std::string> param_names;

        std::ifstream myFile(filename.string());

        if (!myFile.is_open()) throw std::runtime_error("LancFM unable to load file parameter file");
        std::cout << "Processing.." << filename << std::endl;

        // Read the column names
        if (myFile.good())
        {
            // Extract the first line in the file
            std::getline(myFile, header);
            std::stringstream ss(header);


            // Extract each param name
            std::string t;
            int i = 0;
            int columns = 0;
            while (std::getline(ss, t, ','))
            {
                param_names.push_back(t);

                // remember which column special flight condition params are in so we can key on this for insert into data frame
                if (t == "AoA" || t == "Alpha") {
                    alpha_idx = i;
                    printf("Alpha in col %d", i);
                }
                else if (t == "Beta") 
                {
                    beta_idx = i;
                    printf("Beta in col %d", i);
                }
                else if (t == "Mach")
                {
                    mach_idx = i;
                    printf("Mach in col %d", i);
                }
                i++;
            }

            columns = i;

            // Read data, line by line.
            // sometimes VSP .polar does silly stuff and print columns of data without heading, so load only data with headings

            while (std::getline(myFile, line))
            {       
                double row_mach = 0;
                double row_alpha = 0;
                double row_beta = 0;

                std::stringstream ss(line);

                // Keep track of the current column index
                int colIdx = 0;
                std::vector<std::pair<int, double>> val_cache;

                // Extract data into a vector
                while (std::getline(ss, t, ','))
                {
                    if (colIdx == alpha_idx) row_alpha = stod(t);
                    else if (colIdx == beta_idx) row_beta = stod(t);
                    else if (colIdx == mach_idx) row_mach = stod(t);
                    else if (colIdx < columns)
                    {
                        val_cache.push_back(std::make_pair(colIdx, stod(t)));
                    }   // store values in vector
                    colIdx++;
                }

                for (auto& pair : val_cache)
                {
                    insertPolar(param_names[pair.first], row_mach, row_alpha, row_beta, pair.second);
                }

                num_lines++;
            }

            if (!silent) printf("Loaded %d lines", num_lines);
            if(verbose) FMDataLoader::printPolars();
            myFile.close();
            if (!silent) printf("\nFM Loading Complete!");
        }
    }

    void                FMDataLoader::loadVSPcsv(std::filesystem::path filename)
    {
        // load data directly from OpenVSP history.csv
        // sets up elements and polars if found
        // recent versions of VSP appear to write the airframe polars as a separate .polar file, so use FMDataLoader::loadcsv

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

            // deflection (if present) is read from filename
            // to set deflection filename should contain anywhere in its name the following (no quotes) "dElevator&-15.2&"
            // this tells us that the Elevator element is deflected -15.2 degrees in that file
            // 
            // default is 0.0 deflection
            //
            // TODO multiple deflected elements per file?

            double deflected_deg = 0.0;
            std::string deflected_elem;

            // see if the filename identifies a deflected element
            int startpos = filename.filename().string().find("d");
            printf("\nDeflection found at %d", startpos);
            if (startpos != std::string::npos)
            {
                std::string def_str = filename.filename().string().substr(startpos+1);
                std::stringstream def(def_str);

                std::string def_amt_str;
                std::getline(def, deflected_elem, '&');
                printf("\nDeflection elem %s", deflected_elem.c_str());

                std::getline(def, def_amt_str, '&');
                printf("\nDeflection amt %f", deflected_deg);

                deflected_deg = stod(def_amt_str);

                if (verbose) printf("/nDeflection of %f for %s", deflected_deg, deflected_elem.c_str());
            }




            // Extract first column in the file so we can paginate it
            std::string line;

            int linecount = 0;
          //  int total_lines = std::count(std::istreambuf_iterator<char>(myFile), std::istreambuf_iterator<char>(), '\n');  // so we know when we're at the end of the file

            //if (verbose) printf("Total lines = %d", total_lines);

            std::string p = "";      // we'll store each page here temporarily
            std::string cur_page_name = "VSPAERO_History";  // what type of VSP output is this page (elements? totals? other?)

            while (getline(myFile, line))       // first we'll break the text file into consumable Pages delimited by the VSP "Results_Name" line.
            {
                std::string r, t;               // temp string holders
                linecount++;

                std::stringstream ss(line); 
                std::getline(ss, r, ',');       // get first comma-separate token in line
              
                if (r.compare("Results_Name") == 0 ) // store the page start and corresponding VSP page name (e.g. VSPAERO_Comp_Load)
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
            if (p.length() > 0) pages.push_back(make_pair(p, cur_page_name)); // add the last page to the temp page store

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
                double alpha;
                double beta;
                double mach;

                std::vector<double> deflects = { deflected_deg };
                std::vector<double> alphas;
                std::vector<double> betas;
                std::vector<double> machs;

                if (verbose) printf("Processing... \n%s\n", std::get<1>(page).c_str());

                // airframe-total parameter loading
                if (std::get<1>(page).compare("VSPAERO_Polar") == 0) 
                {
                    std::string l;
                    if (verbose) printf("Found a useful page of polars \n%s \n", page.first.c_str());
                    std::vector<std::tuple<std::string, std::vector<double>>>  param_cache;

                    std::stringstream pagestream(page.first);
                    while (getline(pagestream, l, '\n'))  // get the next line from the page
                    {
                        std::string t;
                        if (verbose) printf("\nline \t%s \n", l.c_str());

                        std::stringstream linestream(l);
                        getline(linestream, t, ',');

                        if (verbose) printf("Checking token \t%s\n", t.c_str());
                        if (t == "Alpha")
                        {
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
                            while(getline(linestream, t, ',')) // record betas for this page in order
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
                            while(getline(linestream, t, ',')) // record machs list for this page in order
                            {
                                if (std::find(AFmachs.begin(), AFmachs.end(), stod(t)) == AFmachs.end())
                                {
                                    AFmachs.push_back(stod(t));
                                    if(verbose) printf("adding new mach to polars\t %f\n", stod(t));

                                }
                                machs.push_back(stod(t));
                            }
                        }
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
                            insertPolar(std::get<0>(param), machs[col_indx], alphas[col_indx], betas[col_indx], val);
                            col_indx++;
                        }
                    }
                }


                //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                // element-wise parameter loading..
           
                if (page.second.compare("VSPAERO_Comp_Load")==0)
                {
                    std::string l;
                    std::vector<std::pair<std::string, std::vector<double>>> page_values;
                    std::vector<std::string> elem;  // remember the elements we've found in the results, this should be stable for any single VSP history file.




                    if (verbose) printf("Found a useful result page %s \n", page.first.c_str());

                    // process each line according to what it is
                    std::stringstream pagestream(page.first);

                    while (getline(pagestream, l, '\n'))  // get the next line from the page
                    {
                        std::string t;
                        std::stringstream linestream(l);

                        getline(linestream, t, ',');
                        
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
                                printf("%s\t", x.c_str());
                                elem.push_back(x);  // add this to the vector of param ordering

                                // but is this elem already initialized as an object?
                                bool found = false;
                                for (auto& e : elements)
                                {
                                    if (e.name == x)
                                    {
                                        found = true;
                                    }
                                }
                                if (!found)
                                {
                                    elements.push_back(FMElementData(x)); // initialize the new element
                                    if (verbose) printf("New element identified and inserted %s\n", x.c_str());
                                }
                            }
                        }
                        else if (std::find(paramlist.begin(), paramlist.end(), t) != paramlist.end())  // if I care about the parameter I found, vectorize the line contents
                        {
                            if (verbose) printf("Vectorizing param list %s\n", t.c_str());
                            std::vector<double> vals;
                            std::string x;
                            while (getline(linestream, x, ','))
                            {
                                vals.push_back(stod(x));
                            }

                            page_values.push_back(make_pair(t, vals));   
                        }
                    }

                    // populate page contents into appropriate element objects
                    for (auto & el : elements)
                    {
                        int w = 0;

                        for (auto & e : elem)
                        {
                            if (e == el.name && e != "NONE")
                            {
                                for (auto& v : page_values)
                                {
                                    if (e == deflected_elem)
                                    {
                                        el.insert(v.first, deflected_deg, mach, alpha, beta, v.second[w]);
                                    }
                                    else
                                    {
                                        el.insert(v.first, 0.0, mach, alpha, beta, v.second[w]);
                                    }
                                }
                            }
                            w++;
                        }
                    }
                }
            }
        }



        myFile.close();
        std::sort(AFmachs.begin(), AFmachs.end());
        std::sort(AFalphas.begin(), AFalphas.end());
        std::sort(AFbetas.begin(), AFbetas.end());

        if (!silent)
        {
            printf("FM Loading Complete!\n");
            printf("\nMachs: \n");
            for (auto& x : AFmachs) { printf("%f\t", x); }
            printf("\nAlphas: \n");
            for (auto& x : AFalphas) { printf("%f\t", x); }
            printf("\nBetas: \n");
            for (auto& x : AFbetas) { printf("%f\t", x); }
            printf("\n");

            for (auto& e : elements)
            {
                printf("\n Element: %s\n", e.name.c_str());
                printf("\nDeflections: \n");
                for (auto& x : e.ELdeflects) { printf("%f\t", x); }
                printf("\nMachs: \n");
                for (auto& x : e.ELmachs) { printf("%f\t", x); }
                printf("\nAlphas: \n");
                for (auto& x : e.ELalphas) { printf("%f\t", x); }
                printf("\nBetas: \n");
                for (auto& x : e.ELbetas) { printf("%f\t", x); }
                printf("\n");

                if (verbose) e.printElement();

            }
            if(verbose) for (auto& x : airframe_polars)
            {
                printf("%s, %f, %f, %f, %f = \t%f\n", std::get<0>(x.first).c_str(), std::get<1>(x.first), std::get<2>(x.first), std::get<3>(x.first), x.second);
            }
        }
    }


// Utility functions, accessors, etc.

    FMDataLoader::FMElementData& FMDataLoader::getCompleteElement(std::string element)
    {
        int l = 0;
        int found_index = -1;

        for (auto& m : elements)
        {
            if (m.name.compare(element) == 0) {
                found_index = l;
                break;
            }
            l++;
        }

        if (found_index < 0) { printf("\n %s Not found in data! Be sure to use same element names as the OpenVSP output\n", element.c_str()); }
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

    double              FMDataLoader::getFMParam(std::string element, std::string param, double deflect, double mach, double alpha, double beta) // get list of given parameter (Cl, Cd, etc) values for "element from CSV data
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

    void                FMDataLoader::insertPolar(std::string param, double mach, double alpha, double beta, double value) 
    {
       
        try
        {
            if (!std::get<1>(airframe_polars.insert({ make_tuple(param,  mach, alpha, beta), value })))
            {
                airframe_polars[make_tuple(param, mach, alpha, beta)] = value;
                if (verbose) printf("value was already present, overwritten\n");
            }
            printf("\nInserted new value %f for param %s in airframe polars\n", value, param.c_str());
        }
        catch (const std::exception& e) { 
            if (!silent) printf("%s",e.what());
        }
        // update the list of known dimensions
        if (std::find(AFmachs.begin(), AFmachs.end(), mach) == AFmachs.end()) {         AFmachs.push_back(mach);        std::sort(AFmachs.begin(), AFmachs.end()); }
        if (std::find(AFalphas.begin(), AFalphas.end(), alpha) == AFalphas.end()) {     AFalphas.push_back(alpha);      std::sort(AFalphas.begin(), AFalphas.end()); }
        if (std::find(AFbetas.begin(), AFbetas.end(), beta) == AFbetas.end()) {         AFbetas.push_back(beta);        std::sort(AFbetas.begin(), AFbetas.end()); }

    }

    double              FMDataLoader::getPolar(std::string param, double mach, double alpha, double beta)
    {
        std::vector<double> mach_proxy = getNearest(AFmachs, mach);
        std::vector<double> alph_proxy = getNearest(AFalphas, alpha);
        std::vector<double> beta_proxy = getNearest(AFbetas, beta);

        double a, b;
        double result = 0.0;
        double m_bias = 0.0;
        double a_bias = 0.0;
        double b_bias = 0.0;

        try {
            if (verbose) printf("a test %s %f %f %f %f\n", param.c_str(),  mach_proxy[0], alph_proxy[0], beta_proxy[0]);
            a = airframe_polars.at(make_tuple(param, mach_proxy[0], alph_proxy[0], beta_proxy[0]));

            if (verbose) printf("b test %s %f %f %f %f\n", param.c_str(), mach_proxy[1], alph_proxy[1], beta_proxy[1]);
            b = airframe_polars.at(make_tuple(param, mach_proxy[1], alph_proxy[1], beta_proxy[1]));

            if (a != b)
            {
                if (mach_proxy[0] != mach_proxy[1]) m_bias = (mach_proxy[0] - mach) / (mach_proxy[0] - mach_proxy[1]);
                if (alph_proxy[0] != alph_proxy[1]) a_bias = (alph_proxy[0] - alpha) / (alph_proxy[0] - alph_proxy[1]);
                if (beta_proxy[0] != beta_proxy[1]) b_bias = (beta_proxy[0] - beta) / (beta_proxy[0] - beta_proxy[1]);

                result = a + (b * (m_bias + a_bias + b_bias) / 4.0);
            }
            else
                result = a;

            if (verbose) printf("\nResult %s is %f + (%f * (%f + %f + %f)/4.0) = %f\n", param.c_str(), a, b, m_bias, a_bias, b_bias, result);
        }
        catch (const std::out_of_range& oor) {
            if (!silent) printf("Polar value not found! %s %f %f %f %f\n", param.c_str(), mach, alpha, beta);
        }
        return result;
    }

    std::vector<double> FMDataLoader::getNearest(std::vector<double>& dim, double key)
    {
        std::vector<double>::iterator search_start = dim.begin();
        std::vector<double>::iterator search_end = dim.end();
        std::vector<double>::iterator midpoint;

        std::vector<double>  result;

        if (dim.size() < 1)
        {
            if (verbose) printf("getNearest: input array is empty!\n");
            result = { 0.0, 0.0 };
        }

        if (verbose) printf("getting nearest %f in %d\n", key, dim.size() );

        auto lower = std::lower_bound(dim.begin(), dim.end(), key);
        auto upper = std::upper_bound(dim.begin(), dim.end(), key);

        // boundary protection
        if (lower != dim.begin()) lower = std::prev(lower);
        if (upper == dim.end()) upper = std::prev(upper);
        // if we hit our value exactly (happens often with deflection)
        if (*lower == key) upper = lower;
        if (*upper == key) lower = upper;

        result.push_back(*lower);
        result.push_back(*upper);

        if(verbose) printf("Got nearest %f %f\n", result[0], result[1]);
        return result;
    }

    std::vector<std::string> FMDataLoader::ListElementNames()
    {
        std::vector<std::string> elist;
        for (auto& e : elements)
        {
            elist.push_back(e.name);
        }
        return elist;
    }

    
// element functions

    void                FMDataLoader::FMElementData::insert(std::string param, double deflect, double mach, double alpha, double beta, double value) {

        // Insert 1 new param value into the element's data frame.  If the value already exists, add the new value to the existing.
        // This supports the case where VSPAero breaks up elements into their degen component parts and names them all the same.
        // In that case we expect that the total observable forces and moments of the object are the sum of degen components.


        bool mode = false; // true if inserts to known values should be updates, false if inserts should overwrite prior values
        //printf("Starting new insert into %s\n", this->name);
        if (find(ELdeflects.begin(), ELdeflects.end(), deflect) == ELdeflects.end())  // if new deflection, add it to the unique list of elements
        {
            this->ELdeflects.push_back(deflect);
            std::sort(ELdeflects.begin(), ELdeflects.end());

        }

        if (find(ELmachs.begin(), ELmachs.end(), mach) == ELmachs.end())  // if new mach, add it to the unique list of machs
        {
            this->ELmachs.push_back(mach);
            std::sort(ELmachs.begin(), ELmachs.end());

        }

        if (find(ELalphas.begin(), ELalphas.end(), alpha) == ELalphas.end())  // if new alpha, add it to the unique list of alphas
        {
            this->ELalphas.push_back(alpha);
            std::sort(ELalphas.begin(), ELalphas.end());

        }

        if (find(ELbetas.begin(), ELbetas.end(), beta) == ELbetas.end())  // if new beta, add it to the unique list of betas
        {
            this->ELbetas.push_back(beta);
            std::sort(ELbetas.begin(), ELbetas.end());

        }


        auto curKey = data_frame.find(make_tuple(param, deflect, mach, alpha, beta));

        if (curKey != data_frame.end())  // if this key already exists 
        {
            if (mode)
            {
                //      if (verbose) printf("Updating value %f for param %s in element %s from %f\n", value, param.c_str(), this->name.c_str(), curKey->second);
                curKey->second = curKey->second + value;
            }
            else
            {
                data_frame[make_tuple(param, deflect, mach, alpha, beta)] = value;
            }

        }
        else {
            data_frame[make_tuple(param, deflect, mach, alpha, beta)] = value;
            //   if (verbose) printf("Inserted value %f for param %s in element %s\n", value, param.c_str(), this->name.c_str());
        }
    }

    double              FMDataLoader::FMElementData::lookup(std::string param, double deflect, double mach, double alpha, double beta) {
        // get a single parameter value at a given deflection, mach, alpha, beta
        // par names are same as OpenVSP nomenclature (CL, Cmx, CDi, etc)

        double result = 0.0;
        double a, b;
        double d_bias = 0.0;
        double m_bias = 0.0;
        double a_bias = 0.0;
        double b_bias = 0.0;


        std::vector<double> def_proxy = getNearest(this->ELdeflects, deflect);
        std::vector<double> mach_proxy = getNearest(this->ELmachs, mach);
        std::vector<double> alph_proxy = getNearest(this->ELalphas, alpha);
        std::vector<double> beta_proxy = getNearest(this->ELbetas, beta);

        if (verbose) printf("\nLookup Element %s %s %f %f %f %f\n", this->name.c_str(), param.c_str(), deflect, mach, alpha, beta);
        try {
            a = data_frame.at(std::make_tuple(param, def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]));
            a = data_frame.at(std::make_tuple("CL", 0.0, 0.6, 0.0, 0.0));

        }
        catch (const std::out_of_range& oor) {
            printf("a test failed %s %f %f %f %f\n", param.c_str(), def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]);
            printf(oor.what());
            a = 0.0;
        }
        try {
            b = data_frame.at(make_tuple(param, def_proxy[1], mach_proxy[1], alph_proxy[1], beta_proxy[1]));
        }
        catch (const std::out_of_range& oor) {
            printf("b test failed %s %f %f %f %f\n", param.c_str(), def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]);
            printf(oor.what());
            b = 0.0;
        }

        if (verbose) printf("a test %s %f %f %f %f\n", param.c_str(), def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]);
        if (verbose) printf("b test %s %f %f %f %f\n", param.c_str(), def_proxy[1], mach_proxy[1], alph_proxy[1], beta_proxy[1]);

        if (a != b)
        {
            if (def_proxy[0] != def_proxy[1]) d_bias = (def_proxy[0] - deflect) / (def_proxy[0] - def_proxy[1]);
            if (mach_proxy[0] != mach_proxy[1]) m_bias = (mach_proxy[0] - mach) / (mach_proxy[0] - mach_proxy[1]);
            if (alph_proxy[0] != alph_proxy[1]) a_bias = (alph_proxy[0] - alpha) / (alph_proxy[0] - alph_proxy[1]);
            if (beta_proxy[0] != beta_proxy[1]) b_bias = (beta_proxy[0] - beta) / (beta_proxy[0] - beta_proxy[1]);

            result = a + (b * (d_bias + m_bias + a_bias + b_bias) / 4.0);
        }
        else
        {
            result = a;
            if (verbose) printf("\nOut of mapped range");
        }

        if (verbose) printf("\nElement %s Result %s is %f + (%f * (%f + %f + %f + %f)/4.0) = %f\n", this->name.c_str(), param.c_str(), a, b, d_bias, m_bias, a_bias, b_bias, result);


        /* catch (const std::out_of_range& oor) {
             result = 0.0;
             printf(oor.what());
             if (!silent)
             {
                 printf("\n%s %s value not found! %f %f %f %f\n", this->name.c_str(), param.c_str(), deflect, mach, alpha, beta);
                 printf("a test %s %f %f %f %f\n", param.c_str(), def_proxy[0], mach_proxy[0], alph_proxy[0], beta_proxy[0]);
                 printf("b test %s %f %f %f %f\n", param.c_str(), def_proxy[1], mach_proxy[1], alph_proxy[1], beta_proxy[1]);
             }
         }*/
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

    void                FMDataLoader::FMElementData::printElement() {

        printf("\nElement \t%s\n", this->name.c_str());
        for (auto& p : data_frame)
        {
            printf("%s\t", std::get<0>(p.first).c_str());
            printf("%f\t", std::get<1>(p.first));
            printf("%f\t", std::get<2>(p.first));
            printf("%f\t", std::get<3>(p.first));
            printf("%f\t", std::get<4>(p.first));
            printf("%f\n", p.second);
        }
    }

    void                FMDataLoader::printPolars() {

        for (auto& p : airframe_polars)
        {
            printf("%s\t", std::get<0>(p.first).c_str());
            printf("%f\t", std::get<1>(p.first));
            printf("%f\t", std::get<2>(p.first));
            printf("%f\t", std::get<3>(p.first));
            printf("%f\n", p.second);
        }
    }