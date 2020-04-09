package rmit;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.TreeMap;
import java.util.TreeSet;

import db.Dataset;
import spatialindex.spatialindex.*;
import spatialindex.rtree.*;

public class IKNN {
	private RTree tree;
	private TreeMap<Double, Point> T;
	private HashMap<String, String> trips;
	private Connection conn = null;
	private Dataset ds = null;
	private double lambdaDist;
	public long iotime;
	public long candis;
	
	//TODO
	// Fix lambdaDist / UpperBound / LowerBound
	// Ensure our expected results are then in sorted_candidates
	// Fix main algorithm to ensure we get expected results. 
	
	//With Math.exp(-var) - if var is 0 result is 1. Lower var is higher result is. 

	// compare the distance of two points to given point
	private double compareDistance(Point a, Point b, Point c) {
		double dist_ab = a.getMinimumDistance(b);
		double dist_ac = a.getMinimumDistance(c);

		if (dist_ab > dist_ac) {
			return 1;
		}
		return 0;
	}

	public IKNN(RTree tree, String t, Dataset d, Connection c) throws SQLException {
		// Initializing values
		this.tree = tree;
		this.trips = this.readTrips(t);
		this.ds = d;
		this.conn = c;
	}

	private HashMap<String, String> readTrips(String FilePath) {
		// Hash to store POI and corresponding trajectory IDs
		HashMap<String, String> trips = new HashMap<String, String>();
		try {
			InputStreamReader read = new InputStreamReader(new FileInputStream(FilePath), "utf-8");
			BufferedReader reader = new BufferedReader(read);
			String line;
			String[] arr = null;
			String tmp = null;
			while ((line = reader.readLine()) != null) {
				arr = line.split(",");
				// Parsing value and adding into hash map
				tmp = Double.parseDouble(arr[1]) + "," + Double.parseDouble(arr[2]);

				if (trips.containsKey(tmp)) {
					if (!trips.get(tmp).contains(arr[0])) {
						trips.put(tmp, trips.get(tmp) + "," + arr[0]);
					}
				} else {
					trips.put(tmp, arr[0]);
				}

			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return trips;
	}

	public String computeIKNN(Point[] points, int k) throws SQLException {

		iotime = 0;
		String output = "";
		int lambda = k;

		// stores candidate trajectories (candidate set C)
		HashMap<String, ArrayList<Point>> candidates = new HashMap<String, ArrayList<Point>>();
		// Stores upper bound values (upper bound UBn)
		ArrayList<Point> UB_points = new ArrayList<Point>();
		for (int i = 0; i < points.length; i++) {
			UB_points.add(null);
		}
		int index = 0;
		int iteration = 1;
		int check = 0;
		int counter = 0;
		System.out.println("Number of Query Points:" + points.length);

		// filter the impossible trajectories.
		while (check == 0) {
			System.out.println("Iteration: " + iteration);
			for (int i = 0; i < points.length; i++) {
				// Starts search from first query point and iterates
				index = i;
				ArrayList<Point> S = new ArrayList<Point>();

				long startTime = System.currentTimeMillis();
				// finds lambda nearest points to given query point
				this.getIntersectingPoints(S, lambda, points[i]);
				// System.out.println("Result Point: " + S.size());
				long stopTime = System.currentTimeMillis();
				long elapsedTime = stopTime - startTime;
				iotime += elapsedTime;

				counter = S.size();
				System.out.println("Query Point Index: " + index + " Points: " + S.size());

				for (Point point : S) {

					// 1. stores the farthest point as an upper bound, i.e., update UB_points.
					if (UB_points.get(i) == null) {
						UB_points.set(i, point);
					} else {
						if (point.getMinimumDistance(points[i]) > UB_points.get(i).getMinimumDistance(points[i])) {
							UB_points.set(i, point);
						}
					}

					// 2. finds trajectory ids for each point and stores those trajectories in
					// candidate set, i.e., update candidates.

					String[] tripIDs = trips.get(point.getCoord(0) + "," + point.getCoord(1)).split(",");
					ArrayList<Point> tripTrajectory = null;

					// For each trip ID
					for (int n = 0; n < tripIDs.length; n++) {
						// Load points from database
						tripTrajectory = new ArrayList<Point>();
						tripTrajectory = ds.loadTrajectoryPoints(conn, tripIDs[n]);
						// Add this trajectory to candidates along with tripID
						candidates.put(tripIDs[n], tripTrajectory);
					}
				}
			}

			if (candidates.size() >= k) {

				PriorityQueue<Double> LB = new PriorityQueue<>(Collections.reverseOrder());

				// computes upper bound value using farthest point we found so far for each
				// query point.
				double UB = computeUpperBound(UB_points, points);

				// Compute LB for all trajectories in cadidates
				for (Map.Entry<String, ArrayList<Point>> entry : candidates.entrySet()) {
					// Compute LowerBound
					double tmp_dist = computeLowerBound(entry.getValue(), points);
					LB.add(tmp_dist);
				}

				for (int i = 1; i < k; i++) {
					LB.poll();
					//System.out.println(LB.poll());
				}
				// Choose k-th lower bound value
				double k_LB = LB.peek();

				System.out.println("Current UB: " + UB + " Current k-th LB: " + k_LB);

				if (k_LB >= UB) {
					check = 1;
				}

			}

			lambda += 50;
			iteration++;
			System.out.println("Candidates:" + candidates.size() + "\n");
		}

		System.out.println("End of candidate generation");

		PriorityQueue<Candidate> resultSet = new PriorityQueue<>();// stores top-k results
		PriorityQueue<Candidate> sorted_candidates = new PriorityQueue<>(Collections.reverseOrder());
		// candidates are stored in descending order of upper bound distance.

		for (Map.Entry<String, ArrayList<Point>> entry : candidates.entrySet()) {
			double candidate_ub = computeCandidateUpperBound(entry.getValue(), points, UB_points);
			sorted_candidates.add(new Candidate(candidate_ub, entry.getKey()));
		}
		
		while (sorted_candidates.peek() != null) {
			// 3. scan the candidates and terminate as soon as possible (Algorithm 2).
			for (int i = 1; i <= sorted_candidates.size(); i++) {

				double similarity = computeCandidateDistance(sorted_candidates.peek().getID(), points);

				if (i <= k) {
					Candidate result = sorted_candidates.poll();
					result.setDistance(similarity);
					resultSet.add(result);
				} else {

					if (similarity > resultSet.peek().getDistance()) {
						Candidate result = sorted_candidates.poll();
						result.setDistance(similarity);
						resultSet.poll();
						resultSet.add(result);
					}
					if (i == sorted_candidates.size()
							|| resultSet.peek().getDistance() >= sorted_candidates.peek().getDistance()) {
						break;
					}
				}
			}
			break;
		}

		// format top-k results for printing
		while (resultSet.peek() != null) {
			Candidate can = resultSet.poll();
			output += can.getID() + "\t" + can.getDistance() + "\n";
		}
		this.candis = candidates.size();
		// System.out.println("Iteration: "+iteration+ " Candidates: "+candidates.size()
		// + " Points:" +counter);
		return output;
	}

	// 4 compute the lower bound distance of k-th results
	private double computeLowerBound(ArrayList<Point> p, Point[] points) {
		// find the minDist between each query point and any trajectory point inside
		// lambda distance and sum them for each query point.
		
		//TODO You can't have trajectories which don't have at least one point within lambda distance.
		//That's how we got the fucking trajectories in the first place. 
		//Either lambda distance is wrong or something here is wrong I don't fucking know. 
		
		double dist = 0;
		double minDist = -1;
		double temp;
		for (Point qp : points) {
			for (Point tp : p) {
				temp = tp.getMinimumDistance(qp);
				System.out.println("Actual mindist "+temp);
				if (minDist == -1 && temp < lambdaDist) {
					minDist = temp;
				} else if (temp < minDist && temp < lambdaDist) {
					minDist = temp;
				}
			}
			if (minDist == -1) {
				dist+=0;
			} else {
				dist += Math.exp(-minDist);
				minDist = -1;
			}
			System.out.println("Lamdist: " +this.lambdaDist);
			System.out.println("Mindist: "+ minDist);
			System.out.println("EXPMinDist: " + Math.exp(-minDist) );
		}
		System.out.println("Dist: "+ dist+ "\n");
		return dist;
	}

	// 5 computes upper bound of unseen trajectories using corresponding points
	private double computeUpperBound(ArrayList<Point> p, Point[] points) {
		double dist = 0;
		double temp = 0;
		double lambda = 0;
		for (int i = 0; i < points.length; i++) {
			temp = points[i].getMinimumDistance(p.get(i));
			//System.out.println(temp);
			//System.out.println(Math.exp(-temp));
			dist += Math.exp(-temp);
			lambda+=temp;
		}
		
		this.lambdaDist = (lambda / points.length);
		//System.out.println("Lamdist " + this.lambdaDist);
		return dist;
	}

	// 6 compute the upper bound of candidates
	private double computeCandidateUpperBound(ArrayList<Point> tp, Point[] qp, ArrayList<Point> UB_points) {

		/*
		 * The sum of the distance between points inside lambda and Q points Plus the
		 * remaining points times the UB of unseen traj.
		 */

		double dist = 0;
		double minDist = -1;
		int nonLamCounter = 0;
		for (Point queryPoint : qp) {
			for (Point trajPoint : tp) {
				// Find the minumum distance between a trajectory point and a query point
				if (minDist == -1) {
					minDist = trajPoint.getMinimumDistance(queryPoint);
				} else if (trajPoint.getMinimumDistance(queryPoint) < minDist) {
					minDist = trajPoint.getMinimumDistance(queryPoint);
				}
				//If minimum distance is still outside lambda distance
				if (minDist > this.lambdaDist) {
					minDist = 0;
					nonLamCounter++;
				}
			}
			dist += Math.exp(-minDist);
			minDist = -1;
		}

		dist += (nonLamCounter * this.computeUpperBound(UB_points, qp));

		return dist;
	}

	// 7 Compute actual distance of the trajectory to query points
	public double computeCandidateDistance(String id, Point[] points) throws SQLException {
		double dist = 0;
		double minDist = -1;

		ArrayList<Point> tripTrajectory = new ArrayList<Point>();
		tripTrajectory = ds.loadTrajectoryPoints(conn, id);

		Point testPoint = null;
		
		for (Point queryPoint : points) {
			for (Point trajectoryPoint : tripTrajectory) {
				if (minDist == -1) {
					minDist = trajectoryPoint.getMinimumDistance(queryPoint);
					testPoint = trajectoryPoint;
				} else if (trajectoryPoint.getMinimumDistance(queryPoint) < minDist) {
					minDist = trajectoryPoint.getMinimumDistance(queryPoint);
					testPoint = trajectoryPoint;
				}
			}
			dist += Math.exp(-minDist);
			minDist = -1;
		}
		return dist;
	}

	/*
	 * -- Utility Methods --
	 */

	// Removes the element at the end of a PriorityQueue.
	public PriorityQueue<Candidate> removeLastElement(PriorityQueue<Candidate> p) {
		PriorityQueue<Candidate> newQueue = new PriorityQueue<Candidate>();
		while (p.size() > 1) {
			newQueue.add(p.poll());
		}
		return newQueue;
	}

	// Converts a String array of length 2 to a double array of length 2.
	public double[] StringCoordsToDouble(String[] input) {
		double coord1 = Double.parseDouble(input[0]);
		double coord2 = Double.parseDouble(input[1]);
		double[] coords = { coord1, coord2 };
		return coords;
	}

	public int getIntersectingPoints(ArrayList<Point> s, int lambda, Point p) {
		MyVisitor v = new MyVisitor();
		// finds the nearest lambda points for given query point p
		tree.nearestNeighborQuery(lambda, p, v);

		for (Map.Entry<Integer, IShape> entry : v.answers.entrySet()) {
			IShape value = entry.getValue();
			double[] coord = value.getCenter();
			s.add(new Point(coord));
		}
		return 1;
	}
}
