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
	public long iotime;
	public long candis;

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

				System.out.println("Query Point Index: " + index + " Points: " + S.size());

				for (Point point : S) {
					// Store the farthest point as an upper bound, i.e., update UB_points.
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
						//If this trip is already in candidates
						if (candidates.containsKey(tripIDs[n])) {
							candidates.get(tripIDs[n]).add(point);
						// Add this trajectory to candidates along with tripID
						} else {
							tripTrajectory = new ArrayList<Point>();
							tripTrajectory.add(point);
							candidates.put(tripIDs[n], tripTrajectory);
						}
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
				}
				// Choose k-th lower bound value
				double k_LB = LB.peek();

				System.out.println("Current UB: " + UB + " Current k-th LB: " + k_LB);

				if (k_LB > UB) {
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
					} else {
						sorted_candidates.poll();
					}
					
					//Break condition
					if (i == sorted_candidates.size() || resultSet.peek().getDistance() >= sorted_candidates.peek().getDistance()) {
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
		return output;
	}

	// 4 compute the lower bound distance of k-th results
	private double computeLowerBound(ArrayList<Point> p, Point[] points) {

		double dist = 0;
		double minDist = -1;
		double temp; 
		for (int i = 0; i < points.length; i++) {
			for (int j = 0; j < p.size(); j++) {
				temp = p.get(j).getMinimumDistance(points[i]);
				if (minDist == -1) {
					minDist = temp;
				} else if (temp < minDist) {
					minDist = temp;
				} 
			}
			if (minDist == -1) {
				dist += 0;
			} else {
				dist += Math.exp(-minDist);
				minDist = -1;
			}
		}
		return dist;
	}

	// 5 computes upper bound of unseen trajectories using corresponding points
	private double computeUpperBound(ArrayList<Point> p, Point[] points) {
		double dist = 0;
		double temp = 0;
		for (int i = 0; i < points.length; i++) {
			temp = points[i].getMinimumDistance(p.get(i));
			dist += Math.exp(-temp);
		}
		return dist;
	}

	// 6 compute the upper bound of candidates
	private double computeCandidateUpperBound(ArrayList<Point> tripTrajectory, Point[] points, ArrayList<Point> UB_points) {
		double minDist = -1;
		double dist = 0;
		for (int i = 0; i < points.length; i++) {
			for (Point trajectoryPoint : tripTrajectory) {
				if (minDist == -1) {
					minDist = trajectoryPoint.getMinimumDistance(points[i]);
				} else if (trajectoryPoint.getMinimumDistance(points[i]) < minDist) {
					minDist = trajectoryPoint.getMinimumDistance(points[i]);
				}
			}
			if (minDist > points[i].getMinimumDistance(UB_points.get(i))) {
				minDist = points[i].getMinimumDistance(UB_points.get(i));
			}
			dist += Math.exp(-minDist);
			minDist = -1;
		}
		return dist;
	}

	// 7 Compute actual distance of the trajectory to query points
	public double computeCandidateDistance(String id, Point[] points) throws SQLException {
		double dist = 0;
		double minDist = -1;

		ArrayList<Point> tripTrajectory = new ArrayList<Point>();
		tripTrajectory = ds.loadTrajectoryPoints(conn, id);

		for (Point queryPoint : points) {
			for (Point trajectoryPoint : tripTrajectory) {
				if (minDist == -1) {
					minDist = trajectoryPoint.getMinimumDistance(queryPoint);
				} else if (trajectoryPoint.getMinimumDistance(queryPoint) < minDist) {
					minDist = trajectoryPoint.getMinimumDistance(queryPoint);
				}
			}
			dist += Math.exp(-minDist);
			minDist = -1;
		}
		return dist;
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
