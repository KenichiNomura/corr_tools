#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <cstdlib>
#include <cmath>
#include <cassert>

std::map<std::string,double> CHARGE{ {"O",-2}, {"H",-1} };

struct MDFrame
{
    int natoms;
    float lattice[6];

    std::vector<float> x, y, z;
    std::vector<float> vx, vy, vz;
    std::vector<std::string> name;

	void print()
	{
		std::cout << " # of Atoms : " << natoms << std::endl;
		std::cout << " Lattice Consts : " <<
		lattice[0] << " " << lattice[1] << " " << lattice[2] << " " <<  
		lattice[3] << " " << lattice[4] << " " << lattice[5] << std::endl;

		for(int i=0; i<natoms; i++)
		{
			std::cout << name[i] << " " << x[i] << " " << y[i] << " " << z[i] << " " << 
			vx[i] << " " << vy[i] << " " << vz[i] << std::endl;
		}
	}
};

void print_mdframe(MDFrame &frame)
{
    std::cout << " # of Atoms : " << frame.natoms << std::endl;
    std::cout << " Lattice Consts : " <<
              frame.lattice[0] << " " << frame.lattice[1] << " " << frame.lattice[2] << "\n" <<
              frame.lattice[3] << " " << frame.lattice[4] << " " << frame.lattice[5];

    /*
    	for(int i=0; i<frame.natoms; i++)
    	{
    		std::cout << frame.name[i] << " " << frame.x[i] << " "
    							<< frame.y[i] << " " << frame.z[i] << std::endl;
    	}
    */
}

MDFrame read_single_mdframe(std::ifstream &in)
{
    MDFrame mdframe;
    std::string str;

    std::getline(in,str);
    mdframe.natoms = std::atoi(str.c_str());

    double dummy;
    std::stringstream ss;
    std::getline(in,str);
    ss << str;
    /*
    	ss >> mdframe.lattice[0] >> dummy >> dummy >>
    				dummy >> mdframe.lattice[1] >> dummy >>
    				dummy >> dummy >> mdframe.lattice[2];
    */
    ss >> mdframe.lattice[0] >> mdframe.lattice[1] >> mdframe.lattice[2] >>
	    mdframe.lattice[3] >> mdframe.lattice[4] >> mdframe.lattice[5];

    mdframe.name.resize(mdframe.natoms);
    mdframe.x.resize(mdframe.natoms);
    mdframe.y.resize(mdframe.natoms);
    mdframe.z.resize(mdframe.natoms);
    mdframe.vx.resize(mdframe.natoms);
    mdframe.vy.resize(mdframe.natoms);
    mdframe.vz.resize(mdframe.natoms);

    for (int i=0; i<mdframe.natoms; i++)
    {
        std::string name;
	int id, idummy; 
        float x,y,z,vx,vy,vz, fdummy;

        std::stringstream ss;
        std::getline(in,str);
        ss << str;
        //ss >> name >> x >> y >> z >> vx >> vy >> vz;
        ss >> name >> x >> y >> z >> id >> idummy >> vx >> vy >> vz; 
        //ss >> name >> x >> y >> z >> fdummy >> id >> vx >> vy >> vz; 

        mdframe.name[id-1] = name;
        mdframe.x[id-1] = x;
        mdframe.y[id-1] = y;
        mdframe.z[id-1] = z;
        mdframe.vx[id-1] = vx;
        mdframe.vy[id-1] = vy;
        mdframe.vz[id-1] = vz;
        //std::cout << i << " " << id-1 << " " << name << " " << x << " " << y << " " << z << " " << vx << " " << vy << " " << vz << std::endl;
/*
        mdframe.name.push_back(name);
        mdframe.x.push_back(x);
        mdframe.y.push_back(y);
        mdframe.z.push_back(z);
        mdframe.vx.push_back(vx);
        mdframe.vy.push_back(vy);
        mdframe.vz.push_back(vz);
*/
    }

    return mdframe;
};

struct MSD
{
    std::vector<std::vector<int>> samples;
    std::vector<std::vector<double>> msd;
    std::vector<std::vector<double>> vac, Zt, Gw;
    std::vector<double> vac0;
    std::vector<double> Jt, Jw;
    double J0; 

    std::map<std::string,int> index;

    int total_steps, finalstep_to_add_newref, num_measurements;
    int corr_length, corr_interval;  // [MDsteps]
    double total_steps_fs, corr_length_fs, corr_interval_fs;
    double time_unit_fs; // [fs]

    int num_grids = 1000;
    double freq_unit_fs2eV = 4.13566; // (1fs)^-1 to eV
    double freq_unit_ev2fs = 1.0/freq_unit_fs2eV; // eV to (1fs)^-1
    double max_frequency_ev = 0.05; // [ev]

    MSD(std::vector<MDFrame> & frames,
        double & _corr_length_fs,
        double & _corr_interval_fs,
        double & _time_unit_fs,
	double & _max_frequency_ev)
        :
        corr_length_fs{_corr_length_fs},
        corr_interval_fs{_corr_interval_fs},
        time_unit_fs{_time_unit_fs},
	max_frequency_ev{_max_frequency_ev}
    {
        total_steps = frames.size();
        total_steps_fs = total_steps*time_unit_fs;

        corr_length = (int) (corr_length_fs/time_unit_fs);
        corr_interval = (int) (corr_interval_fs/time_unit_fs);

        assert(corr_length>0);
        assert(corr_interval>0);

        finalstep_to_add_newref = total_steps - corr_length;

        num_measurements = (int) (finalstep_to_add_newref/corr_interval) + 1;
        assert(num_measurements>0);

        auto & frame = frames[0];

        std::cout << "-----------------------------------------------------------------" << std::endl;
        std::cout << " total # of frames : " << frames.size() << std::endl;
        std::cout << " total # of atoms : " << frame.natoms << std::endl;
        std::cout << " lattice parameters : " <<
                  frame.lattice[0] << " " << frame.lattice[1] << " " << frame.lattice[2] << " " <<
                  frame.lattice[3] << " " << frame.lattice[4] << " " << frame.lattice[5] << "\n";

				std::cout << std::endl;

        std::cout << " total steps, corr_length, corr_interval : " 
						<< total_steps << " " << corr_length << " " << corr_interval << " [frames] " << std::endl;
        std::cout << " total steps, corr_length, corr_interval : " 
						<< total_steps_fs << " " << corr_length_fs << " " << corr_interval_fs << " [fs] " << std::endl;
        std::cout << " # of measurements : " << num_measurements << std::endl;

				std::cout << std::endl;

        std::cout << " max_frequency_ev: " << max_frequency_ev*1e3 << " [meV] " << std::endl;
        std::cout << "-----------------------------------------------------------------" << std::endl;

        int index_val = 0;
        for(int i=0; i<frame.natoms; i++)
        {
            auto & name = frame.name[i];
            if(index.find(name) == index.end())
                index[name]=index_val++;
        }

        std::cout << " atom types : " << std::endl;
        for (auto it = index.begin(); it != index.end(); ++it)
            std::cout << "element,id,charge: " << it->first << " " << it->second << " " << CHARGE[it->first] << std::endl ;
            //std::cout << it->first << " => " << it->second << " ";
        std::cout << std::endl;

        msd.resize(index.size());
        vac.resize(index.size());
        vac0.resize(index.size(), 0.0);
        samples.resize(index.size());
        Zt.resize(index.size());
        Gw.resize(index.size());

        Jt.resize(corr_length);
        Jw.resize(num_grids, 0.0);
        J0 = 0.0;

        for(int i=0; i<index.size(); i++)
        {
            msd[i].resize(corr_length, 0.0);
            vac[i].resize(corr_length, 0.0);
            samples[i].resize(corr_length, 0);

            Zt[i].resize(corr_length, 0.0);
            Gw[i].resize(num_grids, 0.0);
        }

    }

    double apply_pbc(double x, double lattice)
    {
        if(x>=0.5*lattice) x-=lattice;
        if(x<-0.5*lattice) x+=lattice;
        return x;
    };

    void get_msd(std::vector<MDFrame> &frames)
    {

        std::vector<MDFrame> ref_frames;

        for(int time=0; time<frames.size(); time++)
        {
            auto & cur_frame = frames[time];

            // add new reference
            if(time%corr_interval==0 && time < finalstep_to_add_newref)
            {
                std::cout << "new reference added. " << time << " / " << frames.size() << std::endl;
                ref_frames.push_back(cur_frame);

		double Jv[3] = {0.0, 0.0, 0.0}; // for current-current correlation

                for(int n=0; n<cur_frame.natoms; n++)
                {
		    auto elem = cur_frame.name[n]; 
                    int type = index[elem];

                    vac0[type] += cur_frame.vx[n]*cur_frame.vx[n] +
                                  cur_frame.vy[n]*cur_frame.vy[n] +
                                  cur_frame.vz[n]*cur_frame.vz[n];

                    Jv[0] = CHARGE[elem]*cur_frame.vx[n];
                    Jv[1] = CHARGE[elem]*cur_frame.vy[n];
                    Jv[2] = CHARGE[elem]*cur_frame.vz[n];
                }

		J0 += Jv[0]*Jv[0] + Jv[1]*Jv[1] + Jv[2]*Jv[2];
            }

            for(int i=0; i<ref_frames.size(); i++)
            {
                auto & ref_frame = ref_frames[i];

                int t = time - i*corr_interval;

                if(t > corr_length) continue;

		double Jv_cur[3] = {0.0, 0.0, 0.0}; // for current-current correlation
                double Jv_ref[3] = {0.0, 0.0, 0.0}; // for current-current correlation

                for(int n=0; n<cur_frame.natoms; n++)
                {
                    double dx = cur_frame.x[n] - ref_frame.x[n];
                    double dy = cur_frame.y[n] - ref_frame.y[n];
                    double dz = cur_frame.z[n] - ref_frame.z[n];
                    //std::cout << dx << " " << dy << " " << dz << std::endl;

                    dx = apply_pbc(dx, cur_frame.lattice[0]);
                    dy = apply_pbc(dy, cur_frame.lattice[1]);
                    dz = apply_pbc(dz, cur_frame.lattice[2]);

                    //std::cout << dx << " " << dy << " " << dz << std::endl;

                    double dr2 = dx*dx + dy*dy + dz*dz;
										//assert(dr2 < 1); // catch too large displacement 

                    double v0vt = cur_frame.vx[n]*ref_frame.vx[n] +
                                  cur_frame.vy[n]*ref_frame.vy[n] +
                                  cur_frame.vz[n]*ref_frame.vz[n];

                    auto elem = cur_frame.name[n];
                    int type = index[cur_frame.name[n]];
				
                    msd[type][t]+=dr2;
                    vac[type][t]+=v0vt;
                    samples[type][t]++;

                    Jv_cur[0] = CHARGE[elem]*cur_frame.vx[n];
                    Jv_cur[1] = CHARGE[elem]*cur_frame.vy[n];
                    Jv_cur[2] = CHARGE[elem]*cur_frame.vz[n];
                    Jv_ref[0] = CHARGE[elem]*ref_frame.vx[n];
                    Jv_ref[1] = CHARGE[elem]*ref_frame.vy[n];
                    Jv_ref[2] = CHARGE[elem]*ref_frame.vz[n];
                }

		Jt[t] += Jv_cur[0]*Jv_ref[0] + Jv_cur[1]*Jv_ref[1] + Jv_cur[2]*Jv_ref[2];
            }
        }

        std::ofstream outfile;
        outfile.open("msd.dat");

        outfile << "Time[fs] ";
        for (auto it = index.begin(); it != index.end(); ++it)
            outfile << "MSD(" << it->first << ") ";
        for (auto it = index.begin(); it != index.end(); ++it)
            outfile << "VAC(" << it->first << ") ";
        for (auto it = index.begin(); it != index.end(); ++it)
            outfile << "Samples(" << it->first << ") ";
        for (auto it = index.begin(); it != index.end(); ++it)
            outfile << "VAC_trunc(" << it->first << ") ";
        outfile << "Jt ";
        outfile << std::endl;

        for(int t=0; t<corr_length; t++)
        {
            for(int type=0; type<samples.size(); type++)
            {
                if(samples[type][t]>0)
                {
                    msd[type][t]/=samples[type][t];
                    vac[type][t]/=vac0[type];
                    Zt[type][t] = vac[type][t]*cos(0.5*M_PI*(double)t/corr_length);
                }
                else
                {
                    msd[type][t]=0.0;
                    vac[type][t]=0.0;
                }
            }

            outfile << t*time_unit_fs << " ";

            for(int type=0; type<samples.size(); type++)
                outfile << msd[type][t] << " ";
            for(int type=0; type<samples.size(); type++)
                outfile << vac[type][t] << " ";
            for(int type=0; type<samples.size(); type++)
                outfile << samples[type][t] << " ";

            for(int type=0; type<samples.size(); type++)
                outfile << Zt[type][t] << " ";

            Jt[t] /= J0;
            outfile << Jt[t];

            outfile << std::endl;

        }

        std::cout << "saved in msd.dat" << std::endl;
        outfile.close();

        outfile.open("dos.dat");

        outfile << "Freq(meV) ";
        for (auto it = index.begin(); it != index.end(); ++it)
            outfile << "DoS(" << it->first << ") ";
        outfile << "DoS(Total) ";
        outfile << "JDoS" << std::endl;

        for(int i=0; i<num_grids; i++) // w
        {
            double w = max_frequency_ev*i/num_grids;

            for(int type=0; type<samples.size(); type++)
            {
                for(int t=0; t<corr_length; t++)
                    Gw[type][i] += std::cos( (2*M_PI*w*freq_unit_ev2fs) * (t*time_unit_fs) )*Zt[type][t]*time_unit_fs;
            }

            for(int t=0; t<corr_length; t++)
                Jw[i] += std::cos( (2*M_PI*w*freq_unit_ev2fs) * (t*time_unit_fs) )*Jt[t]*time_unit_fs;

            double Gw_total = 0.0;
            outfile << w*1e3 << " "; // [meV]
            for(int type=0; type<samples.size(); type++)
            {
                outfile << Gw[type][i] << " ";
                Gw_total += Gw[type][i];
            }
            outfile << Gw_total << " "; 
            outfile << Jw[i] << std::endl;
        }

        std::cout << "saved in dos.dat" << std::endl;
        outfile.close();
    }
};


int main(int argc, char *argv[])
{
		std::vector<MDFrame> mdframes;

		if ( std::string("-h").compare(std::string(argv[1])) == 0)
		{
			std::cout << "help message here" << std::endl;
			exit(0);
		}

		std::string filename = argv[1];
		std::ifstream ifile(filename);

		std::cout << "loading " << filename << std::endl;

		while(!ifile.eof())
			mdframes.push_back(read_single_mdframe(ifile));

		double time_unit, corr_length, corr_interval, max_frequency_ev;

		if (argc == 6)
		{
			corr_length = std::strtod(argv[2], nullptr);
			corr_interval = std::strtod(argv[3], nullptr);
			time_unit = std::strtod(argv[4], nullptr);
				max_frequency_ev = std::strtod(argv[5], nullptr);
		} else {
			std::cout << "use default parameters \n";
			time_unit = 1.0; // 1 [fs]
			corr_length = time_unit*mdframes.size()*0.75; // corr_length is 75% of total data
			corr_interval = corr_length*0.10;  // corr_interval is 25% of correlation length
				max_frequency_ev = 0.1; // 100 [meV]
		}

		MSD msd(mdframes, corr_length, corr_interval, time_unit, max_frequency_ev);
		msd.get_msd(mdframes);
}
