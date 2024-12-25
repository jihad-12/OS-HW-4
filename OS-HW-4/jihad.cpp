#include <iostream>
#include <vector>
#include <algorithm> // for sort
#include <numeric>   // for accumulate
#include <climits>   // for INT_MAX  

using namespace std;

// Define the structure for Process
struct Process {
    int pid;           // Process ID
    int arrival_time;  // Arrival Time
    int cpu_burst;     // CPU Burst Time
    int remaining_time;// Remaining Burst Time
    int start_time;    // Start Time
    int finish_time;   // Finish Time
    int waiting_time;  // Waiting Time
    int turnaround_time; // Turnaround Time
    bool executed;     // Execution status

    // Constructor for Process
    Process(int _pid, int _arrival_time, int _cpu_burst) {
        pid = _pid;
        arrival_time = _arrival_time;
        cpu_burst = _cpu_burst;
        remaining_time = cpu_burst;
        executed = false;
    }
};

// Define the structure for GanttChart
struct GanttChart {
    int pid;           // Process ID
    int start_time;    // Start Time
    int end_time;      // End Time

    // Constructor for GanttChart
    GanttChart(int _pid, int _start_time, int _end_time) {
        pid = _pid;
        start_time = _start_time;
        end_time = _end_time;
    }
};

// Function to print Gantt Chart
void print_gantt_chart(const vector<GanttChart>& gantt_chart) {
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    cout << "| Process | Start Time | End Time   |" << endl;
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
    for (const auto& entry : gantt_chart) {
        cout << "|    " << entry.pid << "    |    " << entry.start_time << "     |    " << entry.end_time << "   |" << endl;
    }
    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~q a1       " << endl;
}

// Function to find waiting time for each process
void findWaitingTime(vector<Process>& processes, vector<int>& wt) {
    int n = processes.size();
    vector<int> rt(n, 0);

    for (int i = 0; i < n; ++i)
        rt[i] = processes[i].cpu_burst;

    int complete = 0, t = 0, minm = INT_MAX;
    int shortest = 0, finish_time;
    bool check = false;

    while (complete != n) {
        for (int j = 0; j < n; ++j) {
            if ((processes[j].arrival_time <= t) && (rt[j] < minm) && rt[j] > 0) {
                minm = rt[j];
                shortest = j;
                check = true;
            } 
        }

        if (!check) {
            ++t;
            continue;
        }

        rt[shortest]--;

        minm = rt[shortest];
        if (minm == 0)
            minm = INT_MAX;

        if (rt[shortest] == 0) {
            ++complete;
            check = false;
            finish_time = t + 1;
            wt[shortest] = finish_time - processes[shortest].cpu_burst - processes[shortest].arrival_time;
            if (wt[shortest] < 0)
                wt[shortest] = 0;
        }
        ++t;
    }
}

// Function to find turnaround time for each process
void findTurnAroundTime(vector<Process>& processes, const vector<int>& wt, vector<int>& tat) {
    int n = processes.size();
    for (int i = 0; i < n; ++i)
        tat[i] = processes[i].cpu_burst + wt[i];
}

// Function to calculate average waiting and turnaround time
void findavgTime(vector<Process>& processes) {
    int n = processes.size();
    vector<int> wt(n, 0), tat(n, 0);

    findWaitingTime(processes, wt);
    findTurnAroundTime(processes, wt, tat);

    cout << " Process\t\t"
        << "BurestTime\t\t"
        << "WaitnigTime\t\t"
        << "TarunAroundTmime\t\t\n";

    for (int i = 0; i < n; ++i) {
        processes[i].waiting_time = wt[i];
        processes[i].turnaround_time = tat[i];
        cout << " " << processes[i].pid << "\t\t"
            << processes[i].cpu_burst << "\t\t " << wt[i]
            << "\t\t " << tat[i] << endl;
    }

    // Calculate average waiting and turnaround time
    float total_wt = accumulate(wt.begin(), wt.end(), 0.0);
    float total_tat = accumulate(tat.begin(), tat.end(), 0.0);

    cout << "\nAverage waiting time = " << total_wt / n << endl;
    cout << "\nAverage turn around time = " << total_tat / n << endl;
}

// Function to find waiting time for each process in (RR)
void findWaitingTimeRR(vector<Process>& processes, vector<int>& wt, int quantum) {
    int n = processes.size();
    vector<int> remaining_time(n, 0);

    for (int i = 0; i < n; ++i)
        remaining_time[i] = processes[i].cpu_burst;

    vector<int> arrival_order;
    for (int i = 0; i < n; ++i)
        arrival_order.push_back(i);

    int t = 0;
    while (!arrival_order.empty()) {
        int pid = arrival_order.front();
        arrival_order.erase(arrival_order.begin());

        if (remaining_time[pid] > 0) {
            if (remaining_time[pid] > quantum) {
                t += quantum;
                remaining_time[pid] -= quantum;
                arrival_order.push_back(pid); // Add process back to the end of the queue
            }
            else {
                t += remaining_time[pid];
                wt[pid] = t - processes[pid].cpu_burst - processes[pid].arrival_time;
                remaining_time[pid] = 0;
            }
        }
    }
}

// Function to find waiting time for each process in (SRTF)
void findWaitingTimeSRTF(vector<Process>& processes, vector<int>& wt) {
    int n = processes.size();
    vector<int> remaining_time(n, 0);

    for (int i = 0; i < n; ++i)
        remaining_time[i] = processes[i].cpu_burst;

    int complete = 0, srtf_current_time = 0;
    vector<int> arrival_order;
    for (int i = 0; i < n; ++i)
        arrival_order.push_back(i);

    while (complete != n) {
        int shortest_process_index = -1;
        int min_remaining_time = INT_MAX;

        for (int i = 0; i < n; ++i) {
            if (processes[i].arrival_time <= srtf_current_time && remaining_time[i] < min_remaining_time && remaining_time[i] > 0) {
                min_remaining_time = remaining_time[i];
                shortest_process_index = i;
            }
        }

        if (shortest_process_index == -1) {
            ++srtf_current_time;
            continue;
        }

        remaining_time[shortest_process_index]--;
        srtf_current_time++;

        if (remaining_time[shortest_process_index] == 0) {
            complete++;
            wt[shortest_process_index] = srtf_current_time - processes[shortest_process_index].arrival_time - processes[shortest_process_index].cpu_burst;
            if (wt[shortest_process_index] < 0)
                wt[shortest_process_index] = 0;
        }
    }
}

int main() {
    // Read processes from input
    int numProcesses;
    cout << "Enter the number of processes: ";
    cin >> numProcesses;

    vector<Process> processes;
    for (int i = 0; i < numProcesses; ++i) {
        int id, arrivalTime, burstTime;
        cout << "Enter details for process " << i + 1 << ":" << endl;
        cout << "PID: ";
        cin >> id;
        cout << "Arrival Time: ";
        cin >> arrivalTime;
        cout << "CPU Burst Time: ";
        cin >> burstTime;
        processes.push_back(Process(id, arrivalTime, burstTime));
    }

    // Read context switch time and quantum for RR from input
    int context_switch, quantum;
    cout << "Enter context switch time: ";
    cin >> context_switch;
    cout << "Enter time quantum for RR: ";
    cin >> quantum;

    // Sort processes by arrival time
    sort(processes.begin(), processes.end(), [](const Process& p1, const Process& p2) {
        return p1.arrival_time < p2.arrival_time;
        });

    // First Come First Serve (FCFS)
    int current_time = 0;
    vector<GanttChart> fcfs_gantt_chart;
    for (auto& process : processes) {
        if (process.arrival_time > current_time) {
            current_time = process.arrival_time;
        }
        process.start_time = current_time;
        current_time += process.cpu_burst; // No context switch needed for FCFS
        process.finish_time = current_time;
        process.waiting_time = process.start_time - process.arrival_time;
        process.turnaround_time = process.finish_time - process.arrival_time;
        fcfs_gantt_chart.push_back(GanttChart(process.pid, process.start_time, process.finish_time));
    }

    // Calculate average waiting and turnaround time for FCFS
    double fcfs_average_waiting_time = accumulate(processes.begin(), processes.end(), 0.0, [](double acc, const Process& process) {
        return acc + process.waiting_time;
        }) / processes.size();
    double fcfs_average_turnaround_time = accumulate(processes.begin(), processes.end(), 0.0, [](double acc, const Process& process) {
        return acc + process.turnaround_time;
        }) / processes.size();
    double fcfs_cpu_utilization = (double)(current_time - processes.front().arrival_time) / current_time;

    // Print results for FCFS
    cout << "Gantt Chart : " << endl;
    cout << "First Come First Serve (FCFS):\n";
    print_gantt_chart(fcfs_gantt_chart);
    cout << "Finish Time\tWaiting Time\tTurnaround Time" << endl;
    for (const auto& process : processes) {
        cout << process.finish_time << "\t\t" << process.waiting_time << "\t\t" << process.turnaround_time << endl;
    }
    cout << "Average Waiting Time = " << fcfs_average_waiting_time << endl;
    cout << "Average Turnaround Time = " << fcfs_average_turnaround_time << endl;
    cout << "CPU Utilization = " << fcfs_cpu_utilization * 100 << "%" << endl;
    cout << endl;

    // Shortest Remaining Time First (SRTF)
    int srtf_current_time = 0;
    vector<GanttChart> srtf_gantt_chart;
    vector<Process> srtf_processes = processes; // Copy of processes for SRTF scheduling

    while (true) {
        int shortest_process_index = -1;
        int min_remaining_time = INT_MAX;

        for (size_t i = 0; i < srtf_processes.size(); ++i) {
            if (srtf_processes[i].remaining_time < min_remaining_time && srtf_processes[i].arrival_time <= srtf_current_time && srtf_processes[i].remaining_time > 0) {
                min_remaining_time = srtf_processes[i].remaining_time;
                shortest_process_index = i;
            }
        }

        if (shortest_process_index == -1) {
            ++srtf_current_time;
            continue;
        }

        srtf_processes[shortest_process_index].remaining_time--;
        srtf_current_time++;

        if (srtf_processes[shortest_process_index].remaining_time == 0) {
            srtf_processes[shortest_process_index].finish_time = srtf_current_time;
        }

        srtf_gantt_chart.push_back(GanttChart(srtf_processes[shortest_process_index].pid, srtf_current_time - 1, srtf_current_time));

        if (srtf_processes[shortest_process_index].remaining_time == 0) {
            ++srtf_current_time;
            srtf_processes[shortest_process_index].turnaround_time = srtf_processes[shortest_process_index].finish_time - srtf_processes[shortest_process_index].arrival_time;
            srtf_processes[shortest_process_index].waiting_time = srtf_processes[shortest_process_index].turnaround_time - srtf_processes[shortest_process_index].cpu_burst;
        }
        if (all_of(srtf_processes.begin(), srtf_processes.end(), [](const Process& p) { return p.remaining_time == 0; }))
            break;
    }

    // Calculate average waiting and turnaround time for SRTF
    vector<int> srtf_waiting_time(srtf_processes.size(), 0);
    findWaitingTimeSRTF(srtf_processes, srtf_waiting_time);

    vector<int> srtf_turnaround_time(srtf_processes.size(), 0);
    findTurnAroundTime(srtf_processes, srtf_waiting_time, srtf_turnaround_time);
    double srtf_average_waiting_time = accumulate(srtf_waiting_time.begin(), srtf_waiting_time.end(), 0.0) / srtf_processes.size();
    double srtf_average_turnaround_time = accumulate(srtf_turnaround_time.begin(), srtf_turnaround_time.end(), 0.0) / srtf_processes.size();
    double srtf_cpu_utilization = (double)accumulate(srtf_processes.begin(), srtf_processes.end(), 0, [](int acc, const Process& process) {
        return acc + process.cpu_burst;
        }) / srtf_current_time;

    // Print results for SRTF
    cout << "Shortest Remaining Time First (SRTF):\n";
    print_gantt_chart(srtf_gantt_chart);
    cout << "Finish Time\tWaiting Time\tTurnaround Time" << endl;
    for (size_t i = 0; i < srtf_processes.size(); ++i) {
        cout << srtf_processes[i].finish_time << "\t\t" << srtf_waiting_time[i] << "\t\t" << srtf_turnaround_time[i] << endl;
    }
    cout << "Average Waiting Time = " << srtf_average_waiting_time << endl;
    cout << "Average Turnaround Time = " << srtf_average_turnaround_time << endl;
    cout << "CPU Utilization = " << srtf_cpu_utilization * 100 << "%" << endl;


    // Round Robin (RR)
    int rr_current_time = 0;
    vector<GanttChart> rr_gantt_chart;
    vector<Process> rr_processes = processes; // Copy of processes for RR scheduling

    while (true) {
        bool all_done = true;
        for (auto& process : rr_processes) {
            if (process.remaining_time > 0) {
                all_done = false;
                int execute_time = min(process.remaining_time, quantum);
                process.remaining_time -= execute_time;
                rr_gantt_chart.push_back(GanttChart(process.pid, rr_current_time, rr_current_time + execute_time));
                rr_current_time += execute_time;
            }
        }
        if (all_done)
            break;
    }

    // Calculate waiting time for RR
    vector<int> rr_waiting_time(rr_processes.size(), 0);
    findWaitingTimeRR(rr_processes, rr_waiting_time, quantum);

    // Calculate average waiting and turnaround time for RR
    vector<int> rr_turnaround_time(rr_processes.size(), 0);
    findTurnAroundTime(rr_processes, rr_waiting_time, rr_turnaround_time);
    double rr_average_waiting_time = accumulate(rr_waiting_time.begin(), rr_waiting_time.end(), 0.0) / rr_processes.size();
    double rr_average_turnaround_time = accumulate(rr_turnaround_time.begin(), rr_turnaround_time.end(), 0.0) / rr_processes.size();
    double rr_cpu_utilization = (double)accumulate(rr_processes.begin(), rr_processes.end(), 0, [](int acc, const Process& process) {
        return acc + process.cpu_burst;
        }) / rr_current_time;

    // Print results for RR
    cout << "Round Robin (RR):\n";
    print_gantt_chart(rr_gantt_chart);
    cout << "Finish Time\tWaiting Time\tTurnaround Time" << endl;
    for (size_t i = 0; i < rr_processes.size(); ++i) {
        cout << rr_processes[i].finish_time << "\t\t" << rr_waiting_time[i] << "\t\t" << rr_turnaround_time[i] << endl;
    }
    cout << "Average Waiting Time = " << rr_average_waiting_time << endl;
    cout << "Average Turnaround Time = " << rr_average_turnaround_time << endl;
    cout << "CPU Utilization = " << rr_cpu_utilization * 100 << "%" << endl;
    cout << endl;

    return 0;
}
