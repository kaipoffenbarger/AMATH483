#ifndef ELEVATOR_HPP
#define ELEVATOR_HPP

#include <iostream>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>
#include <random>
#include <atomic>
#include <vector>
#include <condition_variable>
#include <tuple>
#include <cmath>

using namespace std;

const int NUM_FLOORS = 50;
const int NUM_ELEVATORS = 6;
const int MAX_OCCUPANCY = 10;
const int MAX_WAIT_TIME = 5000; // milliseconds

mutex cout_mtx;
mutex elevator_queue_mtx[NUM_ELEVATORS];
condition_variable elevator_cv[NUM_ELEVATORS];

struct Request {
    int person_id;
    int start_floor;
    int dest_floor;
};

queue<Request> elevator_queues[NUM_ELEVATORS]; // private queue per elevator
vector<int> elevator_positions(NUM_ELEVATORS, 0);
vector<int> elevator_directions(NUM_ELEVATORS, 1); // 1 = up, -1 = down
atomic<int> num_people_serviced(0);
vector<int> global_passengers_serviced(NUM_ELEVATORS, 0);
int npeople;

void safe_print(const string& msg) {
    lock_guard<mutex> lock(cout_mtx);
    cout << msg << endl;
}

int assign_elevator(int start_floor) {
    int min_dist = NUM_FLOORS + 1;
    int best_elevator = 0;
    for (int i = 0; i < NUM_ELEVATORS; ++i) {
        int dist = abs(elevator_positions[i] - start_floor);
        if (dist < min_dist) {
            min_dist = dist;
            best_elevator = i;
        }
    }
    return best_elevator;
}

void person(int id) {
    int start = rand() % NUM_FLOORS;
    int dest = rand() % NUM_FLOORS;
    while (dest == start) {
        dest = rand() % NUM_FLOORS;
    }

    safe_print("Person " + to_string(id) + " wants to go from floor " +
               to_string(start) + " to floor " + to_string(dest));

    Request req{id, start, dest};
    int elevator_id = assign_elevator(start);

    {
        lock_guard<mutex> lock(elevator_queue_mtx[elevator_id]);
        elevator_queues[elevator_id].push(req);
    }

    elevator_cv[elevator_id].notify_one();
}

void elevator(int id) {
    int current_floor = 0;
    int direction = 1;
    vector<Request> passengers;

    while (true) {
        // Drop off
        auto it = passengers.begin();
        while (it != passengers.end()) {
            if (it->dest_floor == current_floor) {
                safe_print("Person " + to_string(it->person_id) +
                           " arrived at floor " + to_string(current_floor));
                ++num_people_serviced;
                ++global_passengers_serviced[id];
                it = passengers.erase(it);
            } else {
                ++it;
            }
        }

        // Pick up requests assigned to this elevator at current floor
        {
            lock_guard<mutex> lock(elevator_queue_mtx[id]);
            queue<Request> remaining;

            while (!elevator_queues[id].empty() && passengers.size() < MAX_OCCUPANCY) {
                Request r = elevator_queues[id].front();
                elevator_queues[id].pop();

                if (r.start_floor == current_floor) {
                    passengers.push_back(r);
                    safe_print("Person " + to_string(r.person_id) +
                               " entered elevator " + to_string(id) +
                               " on floor " + to_string(current_floor));
                } else {
                    remaining.push(r);
                }
            }

            swap(elevator_queues[id], remaining);
        }

        // Exit condition
        {
            lock_guard<mutex> lock(elevator_queue_mtx[id]);
            if (num_people_serviced >= npeople &&
                passengers.empty() &&
                elevator_queues[id].empty()) {
                break;
            }
        }

        // Decide next destination (to next person or destination floor)
        int target_floor = -1;

        {
            lock_guard<mutex> lock(elevator_queue_mtx[id]);
            if (!passengers.empty()) {
                target_floor = passengers.front().dest_floor;
            } else if (!elevator_queues[id].empty()) {
                target_floor = elevator_queues[id].front().start_floor;
            }
        }

        if (target_floor == -1 || target_floor == current_floor) {
            this_thread::sleep_for(chrono::milliseconds(100));
            continue;
        }

        safe_print("Elevator " + to_string(id) +
                   " moving from floor " + to_string(current_floor) +
                   " to floor " + to_string(target_floor));

        int distance = abs(current_floor - target_floor);
        this_thread::sleep_for(chrono::milliseconds(100 * distance));
        current_floor = target_floor;
        elevator_positions[id] = current_floor;
    }

    safe_print("Elevator " + to_string(id) + " has finished servicing all people.");
    safe_print("Elevator " + to_string(id) + " serviced " +
               to_string(global_passengers_serviced[id]) + " passengers.");
}

#endif // ELEVATOR_HPP
