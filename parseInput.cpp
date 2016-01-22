//
// YukawaPlasma class implementation
//   for the parseInput() member function
//
// In-Gee Kim   July, 2015
//

// Private header
#include "YukawaPlasma.h"

int YukawaPlasma::parseInput(ifstream &inputFile)
{
    //
    // YukawaPlasma class member function parseInput()
    //
    // This function parse inputs to assign the member variables,
    //   which defines the YTCP component information, evolution behaviors,
    //   and additional option controls.
    //

    queue<string> componentInput;
    queue<string> evolutionInput;
    queue<string> optionInput;
    queue<string> electronInput;

    int mode;
    string token;

    // provide the initial token
    inputFile >> token;

    while ( token != "&end" ) {
        mode = get_parseMode(token);
        switch ( mode ) {
            case 0: // parse mode for the light component
                inputFile >> token;
                while ( token != "&light" ) {
                    componentInput.push(token);
                    inputFile >> token;
                }
                // operate the component parser
                parseLight(componentInput);
                break;
            case 1: // parse mode for the heavy component
                inputFile >> token;
                while ( token != "&heavy" ) {
                    componentInput.push(token);
                    inputFile >> token;
                }
                // operate the component parser
                parseHeavy(componentInput);
                break;
            case 2: // parse mode for evolution
                inputFile >> token;
                while ( token != "&evolution" ) {
                    evolutionInput.push(token);
                    inputFile >> token;
                }
                // operate the evolution parser
                parseEvolution(evolutionInput);
                break;
            case 3: // parse mode for electron
                inputFile >> token;
                while ( token != "&electron" ) {
                    electronInput.push(token);
                    inputFile >> token;
                }
                // operation the electron parser
                parseElectron(electronInput);
                break;
            case 4: // parse mode for options
                inputFile >> token;
                while ( token != "&options" ) {
                    optionInput.push(token);
                    inputFile >> token;
                }
                // operate the options parser
                parseOptions(optionInput);
                break;
            default: // error in parse mode
                cerr << "Error: YukawaPlasma::parseInput\n"
                    << "\tThe keyword " << token << " does not supported.\n"
                    << "\tCheck the &keyword:\n"
                    << "\t\t&compoent\n"
                    << "\t\t&evolve\n"
                    << "\t\t&electron\n"
                    << "\t\t&options\n"
                    << endl;
                inputFile.close();
                return 10; // error code at the YukawaPlasma constructor
                break;
        } // end of switch ( mode )
        // move to the next token
        inputFile >> token;
    } // end of while ( token != "&end" )

    return 0;
}

int YukawaPlasma::get_parseMode(string &keyword)
{
    int mode;
    if ( keyword == "&light" ) {
        mode = 0;
    } else if ( keyword == "&heavy" ) {
        mode = 1;
    } else if ( keyword == "&evolution" ) {
        mode = 2;
    } else if ( keyword == "&electron" ) {
        mode = 3;
    } else if ( keyword == "&options" ) {
        mode = 4;
    } else {
        // unsupported parse mode
        cerr << "Error: YukawaPlasma::get_ParseMode\n";
        cerr << "\tThe parser mode for the keyword "
            << keyword << " is not supported." << endl;

        exit(10); // error code at the YukawaPlasma constructor
    } // end of if-else

    return mode;
}

void YukawaPlasma::parseLight(queue<string> &list)
{
    // variables for controlling sequence
    const int size = list.size();
    int i = 0;

    // value determinations
    string token;
    string value;

    while (i < size) {
        token = list.front(); // get the first token
        list.pop(); // eliminate the first token from the list
        if ( token == "N_light" ) {
            value = list.front();
            N_light = stoul(value);
            i++;
        } else if ( token == "Gamma_light" ) {
            value = list.front();
            Gamma_light = stod(value);
            i++;
        } else if ( token == "n_light" ) {
            value = list.front();
            n_light = stod(value);
            i++;
        } else if ( token == "m_light" ) {
            value = list.front();
            m_light = stod(value);
            i++;
        } else if ( token == "Z_light" ) {
            value = list.front();
            Z_light = stod(value);
            i++;
        } else { // error
            cerr << "Error: YukawaPlasam::parseLight\n"
                << "\tThe keyword name " << token << " is not supported "
                << "in the &light section.\n"
                << "\t\t N_light\n"
                << "\t\t Gamma_light\n"
                << "\t\t n_light\n"
                << "\t\t m_light\n"
                << "\t\t Z_light\n" << endl;

            exit(10);   // error code at the Yukawa constructor
        } // end of if-else
        list.pop();
        i++;
    } // end of while (i < size)
} // end of YukawaPlasma::parseLight

void YukawaPlasma::parseHeavy(queue<string> &list)
{
    // variables for controlling sequence
    const int size = list.size();
    int i = 0;

    // value determinations
    string token;
    string value;

    while (i < size) {
        token = list.front(); // get the first token
        list.pop(); // eliminate the first token from the list
        if ( token == "N_heavy" ) {
            value = list.front();
            N_heavy = stoul(value);
            i++;
        } else if ( token == "Gamma_heavy" ) {
            value = list.front();
            Gamma_heavy = stod(value);
            i++;
        } else if ( token == "n_heavy" ) {
            value = list.front();
            n_heavy = stod(value);
            i++;
        }else if ( token == "m_heavy" ) {
            value = list.front();
            m_heavy = stod(value);
            i++;
        } else if ( token == "Z_heavy" ) {
            value = list.front();
            Z_heavy = stod(value);
            i++;
        } else { // error
            cerr << "Error: YukawaPlasam::parseHeavy\n"
            << "\tThe keyword name " << token << " is not supported "
            << "in the &heavy section.\n"
            << "\t\t N_heavy\n"
            << "\t\t Gamma_heavy\n"
            << "\t\t n_heavy\n"
            << "\t\t m_heavy\n"
            << "\t\t Z_heavy\n" << endl;

            exit(10);   // error code at the Yukawa constructor
        } // end of if-else
        list.pop();
        i++;
    } // end of while (i < size)
} // end of YukawaPlasma::parseLight

void YukawaPlasma::parseEvolution(queue<string> &list)
{
    // variables for controlling sequence
    const int size = list.size();
    int i = 0;

    // value determinations
    string token;
    string value;

    while (i < size) {
        token = list.front(); // get the first token
        list.pop(); // eliminate the first token from the list
        if ( token == "dt" ) {
            i++;
            value = list.front();
            dt = stod(value);
        } else if ( token == "M_timeSteps" ) {
            i++;
            value = list.front();
            M_timeSteps = stoul(value);
        } else if ( token == "M_snapShots" ) {
            i++;
            value = list.front();
            M_snapShots = stoul(value);
        } else if ( token == "M_eqlb" ) {
            i++;
            value = list.front();
            M_eqlb = stoul(value);
        } else if ( token == "M_preDIH" ) {
            i++;
            value = list.front();
            M_preDIH = stoul(value);
        } else { // error
            cerr << "Error: YukawaPlasma::parseEvolution\n"
                << "\tThe keyword name " << token << " is not supported "
                << "in the &evolution section.\n"
                << "\t\tDelta_t\n"
                << "\t\tM_timeSteps\n"
                << "\t\tM_snapShots\n"
                << "\t\tM_eqlb\n"
                << "\t\tM_preDIH\n" << endl;

            exit(10); // error code at the Yukawa constructor
        } // end of if-else
        list.pop();
        i++;
    } // end of while (i < size)
} // end of YukawaPlasma::parseComponent

void YukawaPlasma::parseElectron(queue<string> &list)
{
    // variables for controlling sequence
    const int size = list.size();
    int i = 0;

    // value determinations
    string token;
    string value;

    while ( i < size ) {
        token = list.front(); // get the first token
        list.pop(); // eliminate the first token from the list
        if ( token == "kappa" )  {
            i++;
            value = list.front();
            kappa = stod(value);
        } else if ( token == "quantumMode" ) {
            i++;
            value = list.front();
            quantumMode = stoul(value);
        } else { // error
            cerr << "Error: YukawaPlasma::parseElectron\n"
            << "\tThe keyword name " << token << " is not supported "
            << "in the &electron section.\n"
            << "\t\tkappa\n"
            << "\t\tquantumMode\n"
            << endl;

            exit(10); // error code at the Yukawa constructor
        } // end of if-else
        list.pop();
        i++;
    } // end of while ( i < size )
} // end of YukawaPlasma::parseOptions

void YukawaPlasma::parseOptions(queue<string> &list)
{
    // variables for controlling sequence
    const int size = list.size();
    int i = 0;

    // value determinations
    string token;
    string value;

    while ( i < size ) {
        token = list.front(); // get the first token
        list.pop(); // eliminate the first token from the list
        if ( token == "initMode" ) {
            i++;
            value = list.front();
            initMode = stoul(value);
        } else if ( token == "finalMode" ) {
            i++;
            value = list.front();
            finalMode = stoul(value);
        } else if ( token == "thermostat" ) {
            i++;
            value = list.front();
            thermostat = stol(value);
        } else { // error
            cerr << "Error: YukawaPlasma::parseOptions\n"
                << "\tThe keyword name " << token << " is not supported "
                << "in the &options section.\n"
                << "\t\tinitMode\n"
                << "\t\tfinalMode\n"
                << "\t\tthermostat\n" << endl;

            exit(10); // error code at the Yukawa constructor
        } // end of if-else
        list.pop();
        i++;
    } // end of while ( i < size )
} // end of YukawaPlasma::parseOptions
