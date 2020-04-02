package rmit;

public class Settings {
	
	public static String project_path = "E:/5.Programming/IDE-Workspaces/Eclipse-General/BDM_A2/";//change to local directory.
	public static String rtree_index_location = project_path + "index/poi";
	public static String trip_ids = project_path + "index/la_trips.txt";
	public static String points_location = project_path + "index/la_points.txt";
	public static String db_name = "triphobo";
	public static String db_user = "root";
	public static String db_pass = "Mike9910!";	
	//sudo apt-get install mysql-server
	//mysql -u root -p
	//create database triphobo
	//exit
	//mysql -u root -p triphobo < /home/sheng/workspace/iknn/src/db/tb_la_dataset.sql
}
