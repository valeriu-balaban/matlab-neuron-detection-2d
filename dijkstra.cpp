#include <iostream>
#include <sstream>
#include <map>
#include <queue>
#include <utility>
#include <math.h>

using namespace std;

using Vector    = pair<float, float>;
using Point     = pair<int16_t, int16_t>;
using Cost      = pair<int16_t, float>;

struct CostPoint{
    Cost    first;
    Point   second;

    CostPoint(Cost c, Point p): first(c), second(p) {};
    
    /*  ATTENTION! Implements > since we need priority queue 
    that returns the lowest value */
    bool operator <(CostPoint const& b) const {
        return !(first < b.first);
    }
};

float dot(Vector const& a, Vector const& b){
    return a.first * b.first + a.second * b.second;
}

bool operator <(Point const& a, Point const& b) {
    if(a.first < b.first)
        return true;

    if(a.first == b.first && a.second < b.second)
        return true;

    return false;
}

bool operator ==(Point const& a, Point const& b) {

    if(a.first == b.first && a.second == b.second)
        return true;
    else
        return false;
}

Point operator +(Point const& a, Point const& b) {
    return Point(a.first + b.first, a.second + b.second);
}

bool operator <(Cost const& a, Cost const& b) {
    if(a.first < b.first)
        return true;

    if(a.first == b.first && a.second < b.second)
        return true;

    return false;
}

/* Checks if the two points are within radius r */
bool withinRadius(Point p, Point q, float r){
    if(abs(p.first - q.first) < r && abs(p.second - q.second) < r)
        return true;
    else
        return false;
}


int main(int argc, char const *argv[]) {
    map<Point, Vector> IMAGE;
    map<Point, float> POINTS;

    uint16_t    MAX_FILLS = 10;
    float       GAMMA = 0.5;

    /* Maximum values for the image coordinates */
    uint32_t R, P, M, N;

    /* Relative coordinates of the 8th neighbors */
    Point neighbors[] = {  
        Point(0, -1), Point(-1, -1), Point(-1, 0), Point(-1,  1),
        Point(0,  1), Point( 1,  1), Point( 1, 0), Point( 1, -1)
    };

    /* Read the size of the data */
    cin >> R >> P >> M >> N >> MAX_FILLS >> GAMMA;

    /*   
        Read the image from the input string as a sequence of lines.
        Each line contains the row and column coordinates of the pixel
        of non zero value, the eigenvalue, and the two components of 
        the eigenvector at that pixel location.
    */
    for(uint32_t i = 0; i < R; i++) {
        Point p;
        Vector v;
        float e;

        cin >> p.first >> p.second >> e >> v.first >> v.second;
        
        /* Set the magnitude of the eigenvector as the eigenvelue */
        v.first *= e;
        v.second *= e;

        /* Save the eigenvector */
        IMAGE[p] = v;
    }

    /*  
        Read the coordinates of the points (neurons) from the input 
        string as a sequence of lines. Each line contains the row and the column
        of the neuron center, followed by the radius of the neuron.
    */
    for(uint32_t i = 0; i < P; i++) {
        
        Point p;
        float r;

        cin >> p.first >> p.second >> r;

        /* Save the location of the neuron */
        POINTS[p] = r;
    }

    /* Print the structures */
    /* cout << "Image:" << endl;
    for(auto x:IMAGE){
        cout << "(" << x.first.first << ", " << x.first.second << ") = ";
        cout << x.second.first << ", " << x.second.second << endl;
    }
    cout << "Neurons:" << endl;
    for(auto x:POINTS){
        cout << "(" << x.first.first << ", " << x.first.second << ") = ";
        cout << x.second << endl;
    } */

    uint16_t p_index = 0;

    /* Trace to path from each point */
    for(auto point:POINTS){
        map<Point, Cost> COST;
        map<Point, Point> PREV;
        map<Point, Point> DEST;
        map<Point, Cost> CONNECTIONS;
        priority_queue<CostPoint> QUEUE;

        Point ps     = point.first;

        cerr << "Processing neuron " << p_index << " of " << P << endl;

        COST[ps] = Cost(0, 0);
        QUEUE.push(CostPoint(Cost(0, 0), ps));

        while(!QUEUE.empty()){
            CostPoint cp = QUEUE.top();
            QUEUE.pop();

            Cost c  = cp.first;
            Point p = cp.second;

            /* All the rest possible paths have larger gaps than the threshold */
            if(c.first > MAX_FILLS)
                break;

            /*  We need to know if we are within the initial neuron radius as
                the cost function inside the initial neuron differs  */
            bool insideOrigin = false;

            bool skipNeighbors = false;

            uint16_t q_index = 0;

            /*  Check if we reached a neuron from the list, and if so, 
                update the connection cost */
            for(auto neighbor:POINTS){
                Point q     = neighbor.first;
                float r     = neighbor.second;

                if(withinRadius(p, q, r)){
                    if(ps == q){
                        insideOrigin = true;
                    } else {
                        /* We reached a neuron */
                        Point link = Point(p_index, q_index);
                        auto it = CONNECTIONS.find(link);

                        skipNeighbors = true;

                        /* Add if not in the connections */
                        if( it == CONNECTIONS.end()){
                            CONNECTIONS[link] = COST[p];
                            DEST[link] = p;

                        /* Update if cost is less */
                        } else if(COST[p] < it->second){
                            it->second = COST[p];
                            DEST[link] = p;
                        }
                    }
                }

                q_index++;
            }

            if(skipNeighbors)
                continue;

            /* Iterate over neightbors */
            for(int i = 0; i < 8; i++){
                Point d         = neighbors[i];
                Point r         = p + d;

                /* Skip if the neighbor is ourtside the image */
                if(r.first < 0 || r.second < 0 || r.first >= M || r.second >= N)
                    continue;

                float norm      = sqrt(dot(d, d));
                Cost  C         = COST[p];
                bool emptyPixel = (IMAGE.find(r) == IMAGE.end());

                /* Avoid inserting elements into the map by querying the values */
                auto itx    = IMAGE.find(p);
                auto ity    = IMAGE.find(r);
                Vector vx   = itx == IMAGE.end() ? Vector (0, 0) : itx->second;
                Vector vy   = ity == IMAGE.end() ? Vector (0, 0) : ity->second;

                Vector dn   = Vector(d.first / norm, d.second / norm);

                /* Empty spaces inside neron do not count towards the limit */
                C.first  += emptyPixel && !insideOrigin ? 1 : 0;
                C.second += GAMMA * (1 - sqrt(dot(vx, vx))); // Eigenvalue contribution
                C.second += (1 - GAMMA) * norm / 2 *
                            (sqrt(1 - abs(dot(dn, vx))) + sqrt(1 - abs(dot(dn, vy)))); // Eigenvector contribution

                /* Update the list of costs */
                auto it = COST.find(r);
                
                if(it == COST.end() || C < it->second){
                    COST[r] = C;
                    PREV[r] = p;
                    QUEUE.push(CostPoint(C, r));
                }
            }
        }

        for(auto x:DEST){
            cout << x.first.first << " " << x.first.second << " ";

            /* Print Path */
            Point t = x.second;

            while(t != ps){
                cout << t.first << " " << t.second << " ";
                t = PREV[t];  
            }
            cout << endl;
        }

        p_index++;
    }
   
    return 0;
}