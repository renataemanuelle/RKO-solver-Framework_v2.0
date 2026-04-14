/************************************************************************************
									IO Functions
*************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <set>
#include <map>
#include <iomanip>
#include <fstream>

#include "Data.h"

/************************************************************************************
 Metodo: WriteSolutionScreen
 Description: Outputs the solution to the screen using the Decoder.
*************************************************************************************/
void WriteSolutionScreen(const char *algorithms[], int numMH, TSol s, 
						 float timeBest, float timeTotal, char instance[], 
						 const TProblemData &data, std::vector <TSol> pool)
{
	printf("\n\n\nRKO: ");
	for (int i=0; i<numMH; i++)
		printf("%s | ", algorithms[i]);
	printf("\nBest MH: %s \nInstance: %s \nsol: ", s.nameMH, instance);
	for (int i=0; i<data.n; i++)
		printf("%.3lf ", s.rk[i]);

	printf("\nofv: %.5lf", s.ofv); 
	printf("\nTotal time: %.3f",timeTotal);
	printf("\nBest time: %.3f\n\n",timeBest);

	// print solution pool 
	printf("\nSolution Pool:\n");
	for (int i = 0; i< (int)pool.size(); i++)
		printf("%.5lf [%s]\n", pool[i].ofv, pool[i].nameMH);

	for (int idx_sel : s.selected_idxs)
	{
    	const Acquisition& a = data.acquisitions[idx_sel];

    	std::cout << "index=" << a.index
              << " | ID=" << a.ID
              << " | sat=" << a.satellite
              << " | time=" << a.time
              << " | score=" << a.score_scenario
              << "\n";
	}	
}

/************************************************************************************
 Metodo: WriteSolution
 Description: Outputs the solution in a txt file using the Decoder.
*************************************************************************************/
void WriteSolution(const char *algorithms[], int numMH, TSol s, 
				   float timeBest, float timeTotal, char instance[], 
				   const TProblemData &data)
{
	char name[256]="../Results/Solutions_RKO";
	strcat(name,".txt");

	// file to write the best solution found
	FILE *solFile;                              

    solFile = fopen(name,"a");

	if (!solFile)
	{
		printf("\n\nFile not found %s!!!",name);
		getchar();
		exit(1);
	}

    fprintf(solFile,"\n\nInstance: %s", instance);
	fprintf(solFile,"\nRKO: ");
	for (int i=0; i<numMH; i++)
		fprintf(solFile,"%s | ", algorithms[i]);

	fprintf(solFile,"\nSol: ");
	for (int i=0; i<data.n; i++)
		fprintf(solFile,"%.3lf ", s.rk[i]);

	fprintf(solFile,"\nofv: %lf", s.ofv);
  	fprintf(solFile,"\nBest time: %.3f",timeBest);
	fprintf(solFile,"\nTotal time:%.3f \n",timeTotal);

	for (int idx_sel : s.selected_idxs)
	{
    	const Acquisition& a = data.acquisitions[idx_sel];

    	std::cout << "index=" << a.index
              << " | ID=" << a.ID
              << " | sat=" << a.satellite
              << " | time=" << a.time
              << " | score=" << a.score_scenario
              << "\n";
	}	

	fclose(solFile);
}

/************************************************************************************
 Metodo: WriteResults
 Description: Outputs the results in a csv file.
*************************************************************************************/
void WriteResults(const char *algorithms[], int numMH, double ofv, 
				  double ofvAverage, std::vector <double> ofvs, float timeBest, 
				  float timeTotal, char instance[])
{
	char name[256]="../Results/Results_RKO";
	strcat(name,".csv");

	FILE *File;
    File = fopen(name,"a");

	if (!File)
	{
		printf("\n\nFile not found %s!!!",name);
		exit(1);
	}

	fprintf(File,"\n%s\t", instance);
	for (int i=0; i<numMH; i++)
		fprintf(File,"%s | ", algorithms[i]);

    fprintf(File,"\t%d", (int)ofvs.size());
    for (unsigned int i=0; i<ofvs.size(); i++){
        fprintf(File,"\t%lf", ofvs[i]);   
	}
	fprintf(File,"\t%lf", ofv);
	fprintf(File,"\t%lf", ofvAverage);
	fprintf(File,"\t%.3f", timeBest);
	fprintf(File,"\t%.3f", timeTotal);

	fclose(File);
}

/************************************************************************************
 Method: EvaluateSolution
 Description: Computes and outputs evaluation metrics for the solution, compatible
              with the EOSPython evaluate() format for cross-framework comparison.
              Prints to screen and writes a CSV file.
*************************************************************************************/
void EvaluateSolution(const TSol &s, const TProblemData &data,
                      float timeBest, float timeTotal, float timeSolver,
                      char instance[], const std::vector<TSol> &pool)
{
    const int n_selected = (int)s.selected_idxs.size();
    if (n_selected == 0)
    {
        printf("\n[Evaluate] No acquisitions selected.\n");
        return;
    }

    // --- Scenario metrics (full instance) ---
    std::set<std::string> all_ids;
    double sum_angle_all = 0, sum_area_all = 0, sum_price_all = 0;
    double sum_sun_all = 0, sum_cloud_all = 0, sum_prio_all = 0;
    for (int i = 0; i < data.n; i++)
    {
        const Acquisition& a = data.acquisitions[i];
        all_ids.insert(a.ID);
        sum_angle_all += a.angle;
        sum_area_all  += a.area;
        sum_price_all += a.price;
        sum_sun_all   += a.sun_elevation;
        sum_cloud_all += a.cloud_cover_real;
        sum_prio_all  += a.priority;
    }
    int total_requests = (int)all_ids.size();
    int total_stereo_pairs = (int)data.stereo_pairs.size();

    // --- Solution metrics (selected only) ---
    std::set<std::string> served_ids;
    double total_score = 0, total_profit = 0, total_area = 0;
    double sum_cloud = 0, sum_angle = 0, sum_prio = 0, sum_sun = 0;
    int cloud_good = 0, cloud_bad = 0;
    int angle_good = 0, angle_bad = 0;
    int prio_counts[5] = {0, 0, 0, 0, 0};
    int stereo_selected = 0, strip_selected = 0;

    for (int idx : s.selected_idxs)
    {
        const Acquisition& a = data.acquisitions[idx];
        served_ids.insert(a.ID);
        total_score  += a.score_scenario;
        total_profit += a.price;
        total_area   += a.area;
        sum_cloud    += a.cloud_cover_real;
        sum_angle    += a.angle;
        sum_prio     += a.priority;
        sum_sun      += a.sun_elevation;

        if (a.cloud_cover_real < 10.0) cloud_good++;
        if (a.cloud_cover_real > 30.0) cloud_bad++;
        if (a.angle <= 10.0) angle_good++;
        if (a.angle >= 30.0) angle_bad++;
        if (a.priority >= 1 && a.priority <= 4) prio_counts[a.priority]++;
        if (a.stereo > 0) stereo_selected++;
        if (a.strips >= 2) strip_selected++;
    }

    int served_requests = (int)served_ids.size();
    double avg_cloud = sum_cloud / n_selected;
    double avg_angle = sum_angle / n_selected;
    double avg_prio  = sum_prio  / n_selected;
    double avg_sun   = sum_sun   / n_selected;

    // --- Stereo pairs completeness ---
    int stereo_pairs_complete = 0;
    for (const auto& [a, b] : data.stereo_pairs)
    {
        if (s.x[a] == 1 && s.x[b] == 1)
            stereo_pairs_complete++;
    }

    // --- Pool diversity ---
    std::set<double> unique_ofvs;
    for (const auto& p : pool) unique_ofvs.insert(p.ofv);
    int pool_diversity = (int)unique_ofvs.size();

    // --- Print to screen ---
    printf("\n");
    printf("================================================================\n");
    printf("  EVALUATION — RKO-EOS\n");
    printf("================================================================\n");

    printf("\n--- SCENARIO (instance) ---\n");
    printf("  %-30s %d\n",    "Requests (unique IDs):",   total_requests);
    printf("  %-30s %d\n",    "Attempts (acquisitions):", data.n);
    printf("  %-30s %d\n",    "Stereo pairs (valid):",    total_stereo_pairs);
    printf("  %-30s %.4f\n",  "Avg angle:",               sum_angle_all / data.n);
    printf("  %-30s %.4f\n",  "Avg area:",                sum_area_all  / data.n);
    printf("  %-30s %.4f\n",  "Avg price:",               sum_price_all / data.n);
    printf("  %-30s %.4f\n",  "Avg sun elevation:",       sum_sun_all   / data.n);
    printf("  %-30s %.4f\n",  "Avg cloud cover (real):",  sum_cloud_all / data.n);
    printf("  %-30s %.4f\n",  "Avg priority:",            sum_prio_all  / data.n);

    printf("\n--- SOLUTION ---\n");
    printf("  %-30s %d\n",    "Acquisitions selected:",   n_selected);
    printf("  %-30s %d\n",    "Unique requests served:",  served_requests);
    printf("  %-30s %.6f\n",  "Total score (scenario):",  total_score);
    printf("  %-30s %.6f\n",  "OFV (negated score):",     s.ofv);
    printf("  %-30s %.4f\n",  "Total profit (price):",    total_profit);
    printf("  %-30s %.4f\n",  "Total area:",              total_area);

    printf("\n  Cloud cover (real):\n");
    printf("    %-28s %.4f\n","Average:",                 avg_cloud);
    printf("    %-28s %d\n",  "Good (< 10):",             cloud_good);
    printf("    %-28s %d\n",  "Bad (> 30):",              cloud_bad);

    printf("\n  Angle:\n");
    printf("    %-28s %.4f\n","Average:",                 avg_angle);
    printf("    %-28s %d\n",  "Good (<= 10):",            angle_good);
    printf("    %-28s %d\n",  "Bad (>= 30):",             angle_bad);

    printf("\n  Priority:\n");
    printf("    %-28s %.4f\n","Average:",                 avg_prio);
    printf("    %-28s %d\n",  "Priority 1:",              prio_counts[1]);
    printf("    %-28s %d\n",  "Priority 2:",              prio_counts[2]);
    printf("    %-28s %d\n",  "Priority 3:",              prio_counts[3]);
    printf("    %-28s %d\n",  "Priority 4:",              prio_counts[4]);

    printf("\n  Sun elevation:\n");
    printf("    %-28s %.4f\n","Average:",                 avg_sun);

    printf("\n--- RKO-SPECIFIC ---\n");
    printf("  %-30s %s\n",    "Best MH:",                 s.nameMH);
    printf("  %-30s %.3f s\n","Time to best:",            timeBest);
    printf("  %-30s %.3f s\n","Solver time (T_solver):",  timeSolver);
    printf("  %-30s %.3f s\n","Total time (T_total):",    timeTotal);
    printf("  %-30s %d / %d\n","Pool diversity:",         pool_diversity, (int)pool.size());
    printf("  %-30s %d\n",    "Stereo pairs complete:",   stereo_pairs_complete);
    printf("  %-30s %d\n",    "Stereo acq selected:",     stereo_selected);
    printf("  %-30s %d\n",    "Strip acq selected:",      strip_selected);

    printf("\n--- SELECTED ACQUISITIONS ---\n");
    printf("  %5s | %6s | %3s | %-19s | %8s | %8s | %8s | %6s | %6s | %4s\n",
           "index", "ID", "sat", "time", "score", "price", "area", "angle", "cloud", "prio");
    printf("  %s\n", std::string(103, '-').c_str());
    for (int idx : s.selected_idxs)
    {
        const Acquisition& a = data.acquisitions[idx];
        printf("  %5d | %6s | %3d | %-19s | %8.4f | %8.1f | %8.1f | %6.2f | %6.2f | %4d\n",
               a.index, a.ID.c_str(), a.satellite, a.time.c_str(),
               a.score_scenario, a.price, a.area, a.angle,
               a.cloud_cover_real, a.priority);
    }
    printf("================================================================\n\n");

    // --- Write CSV file ---
    FILE *csvFile = fopen("../Results/Evaluation_RKO.csv", "w");
    if (!csvFile)
    {
        printf("[Evaluate] Warning: could not write Evaluation_RKO.csv\n");
        return;
    }

    fprintf(csvFile, "section,metric,value\n");

    fprintf(csvFile, "scenario,requests,%d\n",               total_requests);
    fprintf(csvFile, "scenario,attempts,%d\n",               data.n);
    fprintf(csvFile, "scenario,stereo_pairs,%d\n",           total_stereo_pairs);
    fprintf(csvFile, "scenario,avg_angle,%.6f\n",            sum_angle_all / data.n);
    fprintf(csvFile, "scenario,avg_area,%.6f\n",             sum_area_all  / data.n);
    fprintf(csvFile, "scenario,avg_price,%.6f\n",            sum_price_all / data.n);
    fprintf(csvFile, "scenario,avg_sun_elevation,%.6f\n",    sum_sun_all   / data.n);
    fprintf(csvFile, "scenario,avg_cloud_cover,%.6f\n",      sum_cloud_all / data.n);
    fprintf(csvFile, "scenario,avg_priority,%.6f\n",         sum_prio_all  / data.n);

    fprintf(csvFile, "solution,acquisitions,%d\n",           n_selected);
    fprintf(csvFile, "solution,unique_requests_served,%d\n", served_requests);
    fprintf(csvFile, "solution,total_score,%.16f\n",         total_score);
    fprintf(csvFile, "solution,ofv,%.16f\n",                 s.ofv);
    fprintf(csvFile, "solution,total_profit,%.6f\n",         total_profit);
    fprintf(csvFile, "solution,total_area,%.6f\n",           total_area);
    fprintf(csvFile, "solution,avg_cloud_cover,%.6f\n",      avg_cloud);
    fprintf(csvFile, "solution,cloud_good_lt10,%d\n",        cloud_good);
    fprintf(csvFile, "solution,cloud_bad_gt30,%d\n",         cloud_bad);
    fprintf(csvFile, "solution,avg_angle,%.6f\n",            avg_angle);
    fprintf(csvFile, "solution,angle_good_le10,%d\n",        angle_good);
    fprintf(csvFile, "solution,angle_bad_ge30,%d\n",         angle_bad);
    fprintf(csvFile, "solution,avg_priority,%.6f\n",         avg_prio);
    fprintf(csvFile, "solution,priority_1,%d\n",             prio_counts[1]);
    fprintf(csvFile, "solution,priority_2,%d\n",             prio_counts[2]);
    fprintf(csvFile, "solution,priority_3,%d\n",             prio_counts[3]);
    fprintf(csvFile, "solution,priority_4,%d\n",             prio_counts[4]);
    fprintf(csvFile, "solution,avg_sun_elevation,%.6f\n",    avg_sun);

    fprintf(csvFile, "rko,best_mh,%s\n",                     s.nameMH);
    fprintf(csvFile, "rko,time_to_best,%.3f\n",              timeBest);
    fprintf(csvFile, "rko,solver_time,%.3f\n",               timeSolver);
    fprintf(csvFile, "rko,total_time,%.3f\n",                timeTotal);
    fprintf(csvFile, "rko,pool_diversity,%d\n",              pool_diversity);
    fprintf(csvFile, "rko,pool_size,%d\n",                   (int)pool.size());
    fprintf(csvFile, "rko,stereo_pairs_complete,%d\n",       stereo_pairs_complete);
    fprintf(csvFile, "rko,stereo_acq_selected,%d\n",         stereo_selected);
    fprintf(csvFile, "rko,strip_acq_selected,%d\n",          strip_selected);

    fclose(csvFile);
    printf("[Evaluate] Metrics saved to ../Results/Evaluation_RKO.csv\n");
}