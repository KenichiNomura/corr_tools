#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>

class Connection
{
	int num_atoms = 0; 

	typedef std::vector<int> connection_t;

	std::vector<connection_t> connections; 

	public:

	Connection(std::ifstream &fin)
	{
		std::string buffer;
		std::stringstream ss;

		num_atoms=0; 
		while(std::getline(fin, buffer)) {num_atoms++;};
		connections.resize(num_atoms);

		fin.clear();
		fin.seekg(0);
		while(std::getline(fin, buffer))
		{
			ss << buffer;

			int id,iid,jid,type_,num_nbrs; 
			double x_,y_,z_,dr;
			ss >> id >> x_ >> y_ >> z_ >> type_ >> num_nbrs;

			iid = id-1;
			for(int i=0; i<num_nbrs; i++)
			{
				ss >> jid >> dr; 
				connections[iid].push_back(jid);
			}

			//std::cout << buffer << std::endl;; 
			ss.str(""); ss.clear();
		}
	}

	void print()
	{
		for(int i=0; i<num_atoms; i++)
		{
			std::cout << i << " " << connections[i].size() << "   ";
			for(int j=0; j<connections[i].size(); j++)
				std::cout << connections[i][j] << " "; 
			std::cout << std::endl;
		}
		//std::cout << "num_atoms " << num_atoms << std::endl;
	}
};


int main(int argc, char* argv[])
{

	std::string filename(argv[1]);
	std::cout << filename << std::endl;
	std::ifstream fin(filename);

	Connection single_frame(fin); 

	single_frame.print();

	return 0;
}
