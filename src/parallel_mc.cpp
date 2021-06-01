#include "parallel_mc.h"

/*
double ptemp(int num_procs, int id, double action, double alpha, double gamma, double sx, double J, double sflips, IsingLattice2D& lat, bool* switched, std::vector<double>& send_to_meas) {
    //use MPI to perform a parallel tempering step
    //when it comes time to send the lattices, master process will send first
    //and child processes will all receive first.  This might not be the quickest way
    //to do it but it ensures nothing gets locked up waiting to send/receive
    //
    //lattice message tags will be the receiving id

    //return value is the new action value

    //Don't use this function if alpha=0

    //"switched" tells if the lattice is moved or not
    * switched = false;//lattice does not move

    //to begin, action = alpha*c + gamma*gamma_cont + J*J_cont.  we want the variable "action" to be c, so action = (action - gamma*gamma_cont - J*J_cont) / alpha
    MPI_Status Stat[2];
    MPI_Request req[2];
    std::vector<int> lat_buffer_out = lat.get_bool_spins();
    lat_buffer_out.shrink_to_fit();
    std::vector<int> lat_buffer_in(lat.get_N());
    lat_buffer_in.shrink_to_fit();
    double gamma_cont = ((double)lat.get_N()) * (2.0 * sx - 1);
    double J_cont = ((double)lat.get_N()) * (2.0 * sflips - 1);
    double new_action;
    action = (action - gamma * gamma_cont - J * J_cont) / alpha;

    //master process
    //1. receive action from other processes
    //2. test for switches, storing the send-to id in elements of an array
    //3. send the send-to id to each other process
    //4. send lattice if send-to[0] != 0
    if (id == 0) {
        std::vector<double> actions(num_procs);
        std::vector<double> alphas(num_procs);
        std::vector<double> gamma_conts(num_procs);
        std::vector<double> J_conts(num_procs);
        std::vector<double> new_action_list(num_procs);
        std::vector<int> receive_from(num_procs);
        std::vector<int> send_to(num_procs);
        actions[0] = action;
        alphas[0] = alpha;
        gamma_conts[0] = gamma_cont;
        J_conts[0] = J_cont;
        receive_from[0] = 0;
        send_to[0] = 0;
        for (int i = 1; i < num_procs; ++i) {
            MPI_Recv(&(actions[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(alphas[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(gamma_conts[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            MPI_Recv(&(J_conts[i]), 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Stat[0]);
            receive_from[i] = i;
        }

        //test for switches
        //assume for now that action = alpha*C + gamma*gamma_cont
        //start with lowest id and travel up
        for (int i = 0; i < num_procs - 1; ++i) {
            double prob;
            prob = exp((alphas[receive_from[i]] - alphas[receive_from[i + 1]]) * (actions[receive_from[i]] - actions[receive_from[i + 1]]));
            //std::cout << "Action difference: " << actions[receive_from[i]] - actions[receive_from[i + 1]] << "\n";
            //         + 2.0*gamma*(gamma_conts[receive_from[i]] - gamma_conts[receive_from[i + 1]]));this part is probably wrong

            if (drand1_() < prob) {
                //if prob > 1, that means the switching i and i + 1 has a favorable change in action
                //if prob < 1, the switch has an unfavorable change in action but there is still some probability of switching
                receive_from[i + 1] = receive_from[i];
                receive_from[i] = i + 1;
            }
        }
        //invert receive_from to get send_to
        //list new actions for the lattices each process will be receiving
        for (int i = 0; i < num_procs; ++i) {
            send_to[receive_from[i]] = i;
            new_action_list[i] = alphas[i] * actions[receive_from[i]] + gamma * gamma_conts[receive_from[i]] + J * J_conts[receive_from[i]];
        }
        new_action = new_action_list[0];
        //receive_from[i] now gives the id of the lattice that process i will receive
        //message each process and tell them which process they will send their lattice to and which process will receive their lattice
        //tags for send-to id will be the process id, tags for receive-from id will be num_procs + id, tags for new actions will be num_procs + 2*id
        for (int i = 1; i < num_procs; ++i) {
            //send-to
            MPI_Send(&(send_to[i]), 1, MPI_INT, i, i, MPI_COMM_WORLD);
            //receive-from
            MPI_Send(&(receive_from[i]), 1, MPI_INT, i, num_procs + i, MPI_COMM_WORLD);
            //send new action
            MPI_Send(&(new_action_list[i]), 1, MPI_DOUBLE, i, num_procs + 2 * i, MPI_COMM_WORLD);
        }
        //cout << "MPI master process: send_to = " << send_to[0] << ", receive_from = " << receive_from[0] << "\n";

        //send lat_buffer to send_to[0] and receive from receive_from[0]
        if (send_to[0] != 0) {
            MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to[0], send_to[0], MPI_COMM_WORLD, &req[0]);
            MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from[0], 0, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, Stat);
            lat.copy_bool_spins(lat_buffer_in);
            *switched = true; //lattice is moved
        }
        assert(send_to_meas.size() == send_to.size());
        for (int i = 0; i < send_to.size(); ++i) {
            send_to_meas[i] = (double)send_to[i];
        }
    }
    else {
        //child processes
        //1. send action to master process
        //2. receive send-to process and receive-from process
        //3. if send-to[id] != id, send lattice and receive new lattice
        MPI_Send(&action, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&alpha, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&gamma_cont, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);
        MPI_Send(&J_cont, 1, MPI_DOUBLE, 0, id, MPI_COMM_WORLD);


        //cout << "MPI child process sending message action = " << action << ", alpha = " << alpha << "\n";

        int send_to, receive_from;
        MPI_Recv(&send_to, 1, MPI_INT, 0, id, MPI_COMM_WORLD, &Stat[0]);
        MPI_Recv(&receive_from, 1, MPI_INT, 0, num_procs + id, MPI_COMM_WORLD, &Stat[0]);
        MPI_Recv(&new_action, 1, MPI_DOUBLE, 0, num_procs + 2 * id, MPI_COMM_WORLD, &Stat[0]);
        //cout << "MPI child process: send_to = " << send_to << ", receive_from = " << receive_from << "\n";

        //receive lattice and send lattice
        if (send_to != id) {
            MPI_Irecv(&lat_buffer_in.front(), lat_buffer_in.size(), MPI_INT, receive_from, id, MPI_COMM_WORLD, &req[0]);
            //cout << "Process id " << id << " received lattice copy\n";
            MPI_Isend(&lat_buffer_out.front(), lat_buffer_out.size(), MPI_INT, send_to, send_to, MPI_COMM_WORLD, &req[1]);
            MPI_Waitall(2, req, Stat);
            lat.copy_bool_spins(lat_buffer_in);
            *switched = true; //lattice is moved
        }

    }

    return new_action;
}

*/