#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <variant>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <numeric>
#include <math.h>     // sin, cos

using namespace stapl;

const double PI = 3.1415926536;
// Some constants to avoid recalculating them
const double pi2 = 2*PI;      // 2*pi
const double pi4 = 4*PI;      // 4*pi
            

stapl::pVector<std::pair<std::string, std::vector<double>>> read_csv(std::string filename){
    // Reads a CSV file into a vector of <string, vector<float>> pairs where
    // each pair represents <column name, column values>
    using DATA_TYPE = double;
    
    // Create a vector of <string, float vector> pairs to store the result
    stapl::pVector<<std::pair<std::string, std::vector<DATA_TYPE>>> result;

    // Create an input filestream
    std::ifstream myFile(filename);

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname, val_str;
    DATA_TYPE val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);
        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            
            // Initialize and add <colname, float vector> pairs to result
            result.push_back({colname, std::vector<DATA_TYPE> {}});
        }
    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
        // Keep track of the current column index
        int colIdx = 0;
        
        // Extract each float
        while(std::getline(ss, val_str, ',')){
            std::stringstream val_ss(val_str);
            val_ss >> val;
            // Add the current float to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);
            //std::cout << colIdx << " " << val << " " << val_str << " " << result.at(colIdx).second.size() << " line: " << line << std::endl;
            
            // Increment the column index
            colIdx++;
        }  
    }

    // Close file
    myFile.close();
    //std::cout << result.at(1).second.size() << std::endl;
    
    return result;
}

// adapted from https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
void write_csv(std::string filename, std::vector<std::string> header, std::vector<std::string> column_1, std::vector<double> column_2){
    // Make a CSV file with one or more columns of floats
    // Each column of data is represented by the pair <column name, column data>
    //   as std::pair<std::string, std::vector<float>>
    // The data is represented as a vector of these columns
    
    // Create an output filestream object
    std::ofstream myFile(filename);
    
    // Save column names
    for(int j = 0; j < header.size(); ++j)
    {
        myFile << header.at(j);             // every column name
        if(j != header.size() - 1) myFile << ","; // No comma at end of line
    }
    myFile << "\n";
    
    // Save data
    for(int i = 0; i < column_1.size(); ++i)
    {
        myFile << column_1.at(i) << "," << column_1.at(i) << "\n";      // every data vector
    }
    
    // Close the file
    myFile.close();
}

double round_sig(double value, double sig=5){
    /*
    Rounds to a significant number, e.g. to have 5 digits
    :param value: number
    :param sig: int, number of dgits allowed
    :return rounded: float, rounded to the number of sig digits
    */
    if (value == 0.0){ return value; }
    int power = sig-int(floor(log10(abs(value))))-1;
    double multiplier = pow(10,power);
    double rounded = std::ceil(value * multiplier) / multiplier;
    
    return rounded;
}

int bls(stapl::pVector lightcurve) {
    
    namespace fs = std::filesystem;
    std::string path = ".";
    if (argc < 2) {
        std::cout << "No path given, search for waveform files in the current path." << std::endl;
    }
    else {
        path = argv[1];
        std::cout << "Search for waveform files in path " << path << "." << std::endl;
    }
    
    const std::string filename_out = "power_frequencies.csv";
    const std::vector<std::string> result_headers{"Filename", "Frequency [1/s]"};
    std::vector<std::string> power_filenames;
    std::vector<double> power_frequencies;
    
        
    for (const auto & entry : fs::directory_iterator(path)){
        std::string path_filename = entry.path().string();
        
        // split path_filename into path and filename
        std::size_t found = path_filename.find_last_of("/\\");
        std::string path = path_filename.substr(0,found);
        std::string filename = path_filename.substr(found+1);
        
        // are the last 4 characters in lowercase ".csv" and is it not one of the result files?
        if((boost::algorithm::to_lower_copy(path_filename.substr( path_filename.length() - 4 )) == ".csv") and
            (path_filename.find(filename_out) < 0 or path_filename.find(filename_out) > 10000) and
            (filename.find("folded_") != 0) and (filename.find("GLS_") != 0))
        {
            std::cout << "Reading file: " << path_filename << " ... ";
            stapl::pVector<std::pair<std::string, std::vector<double>>> data = read_csv(path_filename);
            //std::cout << data.at(1).second.size() << " " << data.at(0).second.size() << std::endl;
            
            if (data.at(0).second.size() != data.at(1).second.size())
            {
                std::cout << std::endl << "!!! Different number of entries between both rows, please check file. It will be skipped now. ";
                std::cout << "Column " << data.at(0).first << " has " << data.at(0).second.size() << "Entries. ";
                std::cout << "Column " << data.at(1).first << " has " << data.at(1).second.size() << "Entries. " << std::endl;
                continue;
            }
            
            const int number_entries = data.at(1).second.size();
            if (number_entries < 10)        // not enough data
            {
                std::cout << std::endl << "!!! Not enough entries: " << number_entries << " -> will skip." << std::endl;
                continue;
            }
            
            // Check the time unit in the first coloumn
            float time_multi = 1.;    // If not seconds adapt the frequency later
            std::string time_unit = "s";
            if (data.at(0).first != "Seconds") {      // Get units of Time
                if (data.at(0).first == "Minutes") {
                    time_multi = 60.;
                    time_unit = "min";
                }
                else if (data.at(0).first == "Hours") {
                    time_multi = 3600.;
                    time_unit = "h";
                }
                else{
                    std::cout << std::endl << "!!! Unknow time unit: " << data.at(0).first << " , please modify the code. Asuming seconds for now." << std::endl;
                }
            }
            std::cout << "Found entries: " << number_entries;
            
            // put columns into a new vecotr to make code a bit easier to read
            std::vector<double> time_series, time_series_sort, signal_data;
            for (int ii = 0; ii < number_entries; ii++)
            {
                time_series.push_back(data.at(0).second.at(ii));
                time_series_sort.push_back(data.at(0).second.at(ii));
                signal_data.push_back(data.at(1).second.at(ii));
            }
            const double time_range = *std::max_element(time_series.begin(), time_series.end()) - *std::min_element(time_series.begin(), time_series.end());  // returns the numbers
            //int time_range = std::max_element(time_series.begin(), time_series.end()) - std::min_element(time_series.begin(), time_series.end()); // returns the difference of the indexes
            const double min_time_spacing = time_range*1E-9;
            double time_spacing_min = time_range;
            double time_spacing_max = 0.0;
            sort(time_series_sort.begin(), time_series_sort.end());
            for (int ii = 0; ii < number_entries-1; ii++)
            {
                double time_series_diff = time_series_sort[ii+1]-time_series_sort[ii];
                if (time_series_diff > min_time_spacing){
                    // twice the same time in the file would make the diff == 0, which will cause failure at freq_nyqu
                    // using time_range*1E-9 would find 1s signal in 32 years of data length, hence using this instead of "time_series_diff > 0"
                    // using time_range*1E-6 would probably fine too
                    if (time_series_diff < time_spacing_min) {time_spacing_min = time_series_diff;}
                    if (time_series_diff > time_spacing_max) {time_spacing_max = time_series_diff;}
                }
            }
            if (time_spacing_max == 0.0){
                std::cout << std::endl << "!!! All data points have the same time -> will skip." << std::endl;
                continue;
            }
            // Things will still fail if n-1 entries are the same and 1 entry has a different time -> not catched
            
            // Prepare some useful numbers
            const double freq_spacing = 0.2/time_range;       // Frequency spacing, multiplicator could be 0.1 to 0.2
            const double freq_nyqu = 0.5/time_spacing_min;    // Nyquist frquency
            const double freq_max = freq_nyqu/10;             // Useful frequencies: at least 10 data points per period
            
            const double signal_sum = std::accumulate(signal_data.begin(), signal_data.end(), 0.0);
            const double signal_mean = signal_sum / signal_data.size();     // Mean of the signal
            std::vector<double> diff(signal_data.size());
            std::transform(signal_data.begin(), signal_data.end(), diff.begin(), [signal_mean](double x) { return x - signal_mean; });
            const double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            const double signal_variance = sq_sum / diff.size();        // Variance of the signal, do sqrt for standard deviation
            
            std::cout << " , time_range: " << round_sig(time_range,3);
            std::cout << " , time_spacing_min: " << round_sig(time_spacing_min,3) << " , time_spacing_max: " << round_sig(time_spacing_max,3);
            std::cout << " , frequency spacing: " << round_sig(freq_spacing,3);
            std::cout << " , Nyquist frequency: " << round_sig(freq_nyqu,3) << " , maximum frequency to test: " << round_sig(freq_max,4);
            std::cout << " , mean of signal: " << round_sig(signal_mean,3) << " , variance of signal: " << round_sig(signal_variance,3);
            std::cout << " Searching for frequencies ... ";
            
            // Prepare signal - mean(signal)
            std::vector<double> signal_data_mean;
            for (int ii = 0; ii < number_entries-1; ii++)
            {
                signal_data_mean.push_back(signal_data[ii] - signal_mean);
            }
            
            // generalised Lomb-Scargle Periodogram (to not just get the highest frequency)
            std::vector<double> LS_freq, LS_power;
            double freq = freq_spacing;
            while (freq <= freq_max)
            {
                // Some constants to avoid recalculating them
                const double pi4f = pi4 * freq;       // 4*pi*f (=2*omega)
                const double pi2f = pi2 * freq;       // 2*pi*f (=omega)
                // Calculate tau
                double sumsin = 0, sumcos = 0;
                for (int ii = 0; ii < number_entries; ii++)
                {
                    double pi4ft = pi4f*time_series[ii];    // 4*pi*f*t (=2*omega*t)
                    sumsin += sin(pi4ft);
                    sumcos += cos(pi4ft);
                }
                // Eq 34 in https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf
                const double tau = 1/pi4f * atan(sumsin/sumcos);    
                
                // Prepare time-tau
                std::vector<double> time_tau;      // could be done in a lambda as well, but then the vector needs to be copied anyway
                for (int ii = 0; ii < number_entries-1; ii++)
                {
                    time_tau.push_back(time_series[ii] - tau);
                }

                // calculate the generalised LS
                double sumgsin = 0, sumgcos = 0, sumsin2 = 0, sumcos2 = 0;
                for (int ii = 0; ii < number_entries; ii++)
                {
                    const double pi2ft = pi2f*time_tau[ii];    // 4*pi*f*t (=2*omega*t)
                    const double sin_pi2ft = sin(pi2ft);
                    const double cos_pi2ft = cos(pi2ft);
                    sumgsin += signal_data_mean[ii] * sin_pi2ft;
                    sumgcos += signal_data_mean[ii] * cos_pi2ft;
                    sumsin2 += sin_pi2ft * sin_pi2ft;       // sin^2
                    sumcos2 += cos_pi2ft * cos_pi2ft;       // cos^2
                }
                // Eq 33 in https://iopscience.iop.org/article/10.3847/1538-4365/aab766/pdf or Eq 2 in https://arxiv.org/pdf/1502.01344.pdf
                const double PLS = 0.5 * (sumgsin * sumgsin / sumsin2 + sumgcos * sumgcos / sumcos2); 
                
                // Add to results
                LS_freq.push_back(freq);
                LS_power.push_back(PLS);
                
                // Next frequency
                freq += freq_spacing;
            }
            
            /*
            // Put the generalised Lomb-Scargle Periodogram into a file
            std::ofstream myFile(path+"\\GLS_"+filename);
            myFile << "Frequency [1/" << time_unit << "], Power of generalised Lomb-Scargle Periodogram\n";
            for (int ii = 0; ii < LS_power.size(); ii++){
                myFile << LS_freq[ii] << "," << LS_power[ii] << "\n";
            }
            myFile.close();
            */
            
            // Find the maximum power and get the frequency
            auto max_pow_iter = std::max_element(LS_power.begin(), LS_power.end());     // iterator
            auto max_pow_idx = std::distance(LS_power.begin(),max_pow_iter);            // index
            std::cout << "Maximum power in Lomb-Scargle: " << round_sig(LS_power[max_pow_idx],4) << " at frequency " << round_sig(LS_freq[max_pow_idx],4) << " 1/"<< time_unit << std::endl;
            power_frequencies.push_back(round_sig(LS_freq[max_pow_idx]*time_multi,4));   // Frequency in 1/second
            power_filenames.push_back(filename);
            
            /*
            // fold the data with the found frequency
            std::ofstream myFile(path+"\\folded_"+filename);
            myFile << "Phase," << data.at(1).first << "\n";
            for (int ii = 0; ii < signal_data.size(); ii++){
                double time_freq = time_series[ii] * LS_freq[max_pow_idx];
                myFile << time_freq - floor(time_freq+0.5) << "," << signal_data[ii] << "\n";
            }
            myFile.close();
            */
            
        }
    
    }
    
    // Get the indices of the sorted frequencies
    std::vector<int> indices(power_frequencies.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
                [&](int A, int B) -> bool {
                    return power_frequencies[A] < power_frequencies[B];
                });

    std::cout << "Writing result file: " << filename_out << " ... ";
    // Creating result csv file using indices[ii] as index to sort by frequency
    std::ofstream myFile(filename_out);
    myFile << result_headers[0] << ", " << result_headers[1] << "\n";
    for (int ii = 0; ii < indices.size(); ii++){
        myFile << power_filenames[indices[ii]] << "," << power_frequencies[indices[ii]] << "\n";
    }
    myFile.close();
    //write_csv(filename_out, result_headers, waveform_filenames, waveform_frequencies);    // not able to have text and double that way
    std::cout << "Done" << std::endl;
    
    return 0;
}
