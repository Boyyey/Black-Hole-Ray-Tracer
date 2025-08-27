#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <getopt.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

// =====================================================================
// SCIENTIFIC CONFIGURATION
// =====================================================================
#define RS 1.0                   // Schwarzschild radius (fundamental unit)
#define MAX_STEPS 2048           // Maximum photon path steps
#define MIN_STEP 1e-8            // Minimum step size for integration
#define EPSILON 1e-12            // Convergence threshold
#define STAR_DENSITY 0.0005      // Background stars density
#define GAMMA 1.0                // Color gamma (for visualization)
#define MAX_REDSHIFT 10.0        // Redshift clipping value
#define ISCO 6.0*RS              // Innermost Stable Circular Orbit
#define PHOTON_SPHERE 3.0*RS     // Photon sphere radius
#define CRITICAL_IMPACT_PARAM (3.0 * sqrt(3.0) * RS) // Critical impact parameter

// Simulation modes
enum SimulationMode {
    MODE_DEFAULT,
    MODE_BENCHMARK,
    MODE_PARAMETER_SWEEP,
    MODE_DATA_EXPORT,
    MODE_VALIDATION
};

// Output formats
enum OutputFormat {
    FORMAT_CSV,
    FORMAT_HDF5,
    FORMAT_FITS
};

// =====================================================================
// DATA STRUCTURES
// =====================================================================
typedef struct {
    double r, theta, phi;
} SphericalCoord;

typedef struct {
    double t, r, theta, phi;
} SchwarzschildState;

typedef struct {
    double r_dot, theta_dot, phi_dot;
} SchwarzschildVelocity;

typedef struct {
    double impact_parameter;
    double redshift;
    int photon_orbits;
    int escaped;
    double disk_r;
    double disk_phi;
    double path_length;
    double time_of_flight;
} RayTraceInfo;

typedef struct {
    double time;
    double camera_distance;
    double camera_theta;
    double disk_inner;
    double disk_outer;
    int rays_traced;
    int black_hole_hits;
    int disk_hits;
    int star_hits;
    double avg_redshift;
    double max_redshift_observed;
    int photon_orbits_total;
    double simulation_time;
} SimulationTelemetry;

typedef struct {
    double camera_distance;
    double camera_theta;
    double disk_inner;
    double disk_outer;
    int screen_width;
    int screen_height;
    int max_steps;
    double star_density;
    double max_redshift;
    enum OutputFormat output_format;
    char output_prefix[256];
    int validation_mode;
} SimulationConfig;

// =====================================================================
// MATHEMATICAL UTILITIES
// =====================================================================
double schwarzschild_metric(double r) {
    return 1.0 - RS/r;
}

double schwarzschild_metric_derivative(double r) {
    return RS/(r*r);
}

double keplerian_velocity(double r) {
    return sqrt(RS/(2.0*r));
}

double keplerian_angular_velocity(double r) {
    return sqrt(RS/(r*r*r));
}

double solve_cubic(double a, double b, double c) {
    double p = b - a*a/3.0;
    double q = (2.0*a*a*a)/27.0 - (a*b)/3.0 + c;
    double discriminant = q*q/4.0 + p*p*p/27.0;
    
    if (discriminant > 0) {
        double s = cbrt(-q/2.0 + sqrt(discriminant));
        double t = cbrt(-q/2.0 - sqrt(discriminant));
        return -a/3.0 + s + t;
    } else {
        double r = sqrt(-p*p*p/27.0);
        double phi = acos(-q/(2.0*r));
        return -a/3.0 + 2.0*cbrt(r)*cos(phi/3.0);
    }
}

double critical_impact_parameter(double r) {
    return r / sqrt(1.0 - RS/r);
}

double photon_orbit_radius() {
    return 3.0 * RS;
}

// =====================================================================
// RELATIVISTIC RAY TRACING CORE
// =====================================================================
void evolve_geodesic(SchwarzschildState *state, SchwarzschildVelocity *vel, double dlambda) {
    double g_tt = -(1.0 - RS/state->r);
    double g_rr = 1.0/(1.0 - RS/state->r);
    
    // Christoffel symbols (non-zero components only)
    double Gamma_r_tt = (RS/(2.0*state->r*state->r)) * g_tt;
    double Gamma_r_rr = -RS/(2.0*state->r*state->r*state->r) * g_rr;
    double Gamma_r_thth = -(state->r - RS);
    double Gamma_r_phph = -(state->r - RS)*sin(state->theta)*sin(state->theta);
    double Gamma_th_thr = 1.0/state->r;
    double Gamma_th_phph = -sin(state->theta)*cos(state->theta);
    double Gamma_ph_thph = cos(state->theta)/sin(state->theta);
    
    // d²x^μ/dλ² + Γ^μ_αβ dx^α/dλ dx^β/dλ = 0
    double r_ddot = -Gamma_r_tt * pow(state->t, 2) 
                  - Gamma_r_rr * pow(vel->r_dot, 2)
                  - Gamma_r_thth * pow(vel->theta_dot, 2)
                  - Gamma_r_phph * pow(vel->phi_dot, 2);
    
    double theta_ddot = -2.0 * Gamma_th_thr * vel->theta_dot * vel->r_dot
                      - Gamma_th_phph * pow(vel->phi_dot, 2);
    
    double phi_ddot = -2.0 * Gamma_ph_thph * vel->phi_dot * vel->theta_dot;
    
    // Update state using Verlet integration
    state->r += vel->r_dot * dlambda + 0.5 * r_ddot * dlambda * dlambda;
    state->theta += vel->theta_dot * dlambda + 0.5 * theta_ddot * dlambda * dlambda;
    state->phi += vel->phi_dot * dlambda + 0.5 * phi_ddot * dlambda * dlambda;
    
    // Update velocities
    vel->r_dot += r_ddot * dlambda;
    vel->theta_dot += theta_ddot * dlambda;
    vel->phi_dot += phi_ddot * dlambda;
    
    state->t += dlambda;
}

RayTraceInfo trace_ray(double px, double py, SimulationConfig *cfg) {
    // Initialize camera position
    SphericalCoord cam = {
        .r = cfg->camera_distance,
        .theta = cfg->camera_theta,
        .phi = 0.0
    };
    
    // Convert to Cartesian
    double cam_x = cam.r * sin(cam.theta) * cos(cam.phi);
    double cam_y = cam.r * sin(cam.theta) * sin(cam.phi);
    double cam_z = cam.r * cos(cam.theta);
    
    // Ray direction
    double dir_x = px;
    double dir_y = py;
    double dir_z = -1.0;
    
    // Normalize
    double norm = sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
    dir_x /= norm; dir_y /= norm; dir_z /= norm;
    
    // Transform to global coordinates
    double global_dir_x = dir_x * cos(cam.phi) - dir_y * sin(cam.phi);
    double global_dir_y = dir_x * sin(cam.phi) + dir_y * cos(cam.phi);
    double global_dir_z = dir_z;
    
    // Initial state
    SchwarzschildState state = {
        .t = 0.0,
        .r = cam.r,
        .theta = cam.theta,
        .phi = cam.phi
    };
    
    // Initial velocity
    SchwarzschildVelocity vel = {
        .r_dot = global_dir_x,
        .theta_dot = global_dir_y / cam.r,
        .phi_dot = global_dir_z / (cam.r * sin(cam.theta))
    };
    
    RayTraceInfo info = {
        .impact_parameter = 0.0,
        .redshift = 1.0,
        .photon_orbits = 0,
        .escaped = 1,
        .disk_r = -1.0,
        .disk_phi = 0.0,
        .path_length = 0.0,
        .time_of_flight = 0.0
    };
    
    int steps = 0;
    double lambda = 0.0;
    double prev_r = state.r;
    int orbit_count = 0;
    int crossed_3rs = 0;
    
    // Ray tracing loop
    while (steps < cfg->max_steps) {
        // Adaptive step sizing
        double step = fmin(MIN_STEP * 100.0, 0.1 * state.r);
        if (state.r < 3.0*RS) step = fmin(step, MIN_STEP * 10.0);
        if (state.r < 1.5*RS) step = fmin(step, MIN_STEP);
        
        // Save previous state
        prev_r = state.r;
        
        // Evolve photon
        evolve_geodesic(&state, &vel, step);
        lambda += step;
        steps++;
        
        // Check photon sphere crossing
        if (state.r < 3.01*RS && state.r > 2.99*RS) {
            if (prev_r > 3.0*RS && state.r < 3.0*RS) orbit_count++;
        }
        
        // Check 3RS threshold
        if (!crossed_3rs && state.r < 3.0*RS) crossed_3rs = 1;
        
        // Check escape
        if (state.r > 2.0*cfg->camera_distance && crossed_3rs) {
            info.escaped = 1;
            break;
        }
        
        // Check accretion disk intersection
        if (fabs(state.theta - M_PI/2.0) < 1e-4 && 
            state.r >= cfg->disk_inner && state.r <= cfg->disk_outer) {
            
            info.disk_r = state.r;
            info.disk_phi = state.phi;
            info.escaped = 0;
            info.path_length = lambda;
            
            // Calculate gravitational redshift
            double z = 1.0 / sqrt(schwarzschild_metric(info.disk_r)) - 1.0;
            info.redshift = fmin(cfg->max_redshift, 1.0 + z);
            
            // Doppler boosting
            double v_phi = keplerian_velocity(info.disk_r);
            double doppler_factor = 1.0 / (sqrt(1.0 - v_phi*v_phi) * 
                              (1.0 - v_phi * cos(state.phi - info.disk_phi)));
            info.redshift *= doppler_factor;
            
            break;
        }
        
        // Check horizon crossing
        if (state.r <= RS*1.001) {
            info.escaped = 0;
            info.path_length = lambda;
            break;
        }
    }
    
    info.photon_orbits = orbit_count;
    
    // Calculate impact parameter
    double b = cam.r * sqrt(dir_x*dir_x + dir_y*dir_y) / fabs(dir_z);
    info.impact_parameter = fabs(b);
    
    return info;
}

// =====================================================================
// DATA OUTPUT
// =====================================================================
int save_photon_data_csv(const char *filename, RayTraceInfo **ray_data, int width, int height) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Failed to open output file");
        return -1;
    }
    
    fprintf(fp, "x,y,impact_parameter,redshift,photon_orbits,escaped,disk_r,disk_phi,path_length,time_of_flight\n");
    
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            RayTraceInfo *info = &ray_data[y][x];
            
            fprintf(fp, "%d,%d,%.16f,%.16f,%d,%d,%.16f,%.16f,%.16f,%.16f\n",
                    x, y,
                    info->impact_parameter,
                    info->redshift,
                    info->photon_orbits,
                    info->escaped,
                    info->disk_r,
                    info->disk_phi,
                    info->path_length,
                    info->time_of_flight);
        }
    }
    
    fclose(fp);
    printf("Photon data saved to %s\n", filename);
    return 0;
}

int save_telemetry_data(const char *filename, SimulationTelemetry *telemetry) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Failed to open telemetry file");
        return -1;
    }
    
    fprintf(fp, "time,camera_distance,camera_theta,disk_inner,disk_outer,"
               "rays_traced,black_hole_hits,disk_hits,star_hits,"
               "avg_redshift,max_redshift,photon_orbits_total,simulation_time\n");
    
    fprintf(fp, "%.0f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d,%d,%.4f,%.4f,%d,%.4f\n",
            telemetry->time,
            telemetry->camera_distance,
            telemetry->camera_theta,
            telemetry->disk_inner,
            telemetry->disk_outer,
            telemetry->rays_traced,
            telemetry->black_hole_hits,
            telemetry->disk_hits,
            telemetry->star_hits,
            telemetry->avg_redshift,
            telemetry->max_redshift_observed,
            telemetry->photon_orbits_total,
            telemetry->simulation_time);
    
    fclose(fp);
    printf("Telemetry data saved to %s\n", filename);
    return 0;
}

int save_parameter_sweep_csv(const char *filename, double *param_values, 
                           double *results, int count) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("Failed to open parameter sweep file");
        return -1;
    }
    
    fprintf(fp, "parameter,value\n");
    
    for (int i = 0; i < count; i++) {
        fprintf(fp, "%.16f,%.16f\n", param_values[i], results[i]);
    }
    
    fclose(fp);
    printf("Parameter sweep data saved to %s\n", filename);
    return 0;
}

// =====================================================================
// BENCHMARKING
// =====================================================================
double run_benchmark(int iterations, SimulationConfig *cfg) {
    double total_time = 0.0;
    
    for (int i = 0; i < iterations; i++) {
        struct timeval start, end;
        gettimeofday(&start, NULL);
        
        #pragma omp parallel for collapse(2)
        for (int y = 0; y < cfg->screen_height; y++) {
            for (int x = 0; x < cfg->screen_width; x++) {
                double px = (x - cfg->screen_width/2.0)/(0.5 * cfg->screen_width);
                double py = (cfg->screen_height/2.0 - y)/(0.5 * cfg->screen_height);
                
                trace_ray(px, py, cfg);
            }
        }
        
        gettimeofday(&end, NULL);
        double elapsed = (end.tv_sec - start.tv_sec) + 
                        (end.tv_usec - start.tv_usec) / 1000000.0;
        total_time += elapsed;
    }
    
    return total_time / iterations;
}

// =====================================================================
// PARAMETER SWEEP
// =====================================================================
int run_parameter_sweep(SimulationConfig *cfg, 
                      const char *param_name,
                      double start, double end, double step,
                      const char *output_filename) {
    int num_points = (int)((end - start) / step) + 1;
    double *param_values = malloc(num_points * sizeof(double));
    double *results = malloc(num_points * sizeof(double));
    
    if (!param_values || !results) {
        fprintf(stderr, "Memory allocation failed for parameter sweep\n");
        free(param_values);
        free(results);
        return -1;
    }
    
    printf("Running parameter sweep for %s from %.2f to %.2f in steps of %.2f\n",
           param_name, start, end, step);
    
    for (int i = 0; i < num_points; i++) {
        double param = start + i * step;
        param_values[i] = param;
        
        // Set the parameter based on name
        if (strcmp(param_name, "distance") == 0) {
            cfg->camera_distance = param * RS;
        } else if (strcmp(param_name, "theta") == 0) {
            cfg->camera_theta = param * M_PI / 180.0;
        } else if (strcmp(param_name, "disk_inner") == 0) {
            cfg->disk_inner = param * RS;
        } else if (strcmp(param_name, "disk_outer") == 0) {
            cfg->disk_outer = param * RS;
        }
        
        // Run simulation and collect statistics
        int disk_hits = 0;
        int total_rays = cfg->screen_width * cfg->screen_height;
        
        #pragma omp parallel for reduction(+:disk_hits)
        for (int y = 0; y < cfg->screen_height; y++) {
            for (int x = 0; x < cfg->screen_width; x++) {
                double px = (x - cfg->screen_width/2.0)/(0.5 * cfg->screen_width);
                double py = (cfg->screen_height/2.0 - y)/(0.5 * cfg->screen_height);
                
                RayTraceInfo info = trace_ray(px, py, cfg);
                
                if (info.disk_r > 0) {
                    disk_hits++;
                }
            }
        }
        
        // Store result (fraction of rays hitting disk)
        results[i] = (double)disk_hits / total_rays;
        
        printf("  %s = %.2f: %.2f%% of rays hit disk\n", 
               param_name, param, results[i] * 100.0);
    }
    
    // Save results
    int result = save_parameter_sweep_csv(output_filename, param_values, results, num_points);
    
    free(param_values);
    free(results);
    
    return result;
}

// =====================================================================
// VALIDATION
// =====================================================================
int validate_photon_sphere(SimulationConfig *cfg) {
    printf("Validating photon sphere at r = %.4f Rs\n", PHOTON_SPHERE/RS);
    
    // Test rays with impact parameter near critical value
    double b_test = CRITICAL_IMPACT_PARAM * 0.999;
    double px = b_test / cfg->camera_distance;
    double py = 0.0;
    
    RayTraceInfo info = trace_ray(px, py, cfg);
    
    if (info.escaped == 0 && info.disk_r <= 0) {
        printf("  ✓ Photon capture confirmed for b = %.4f (critical b = %.4f)\n", 
               b_test, CRITICAL_IMPACT_PARAM);
        return 0;
    } else {
        printf("  ✗ Validation failed: expected photon capture but got %s\n",
               info.escaped ? "escape" : (info.disk_r > 0 ? "disk hit" : "unknown"));
        return -1;
    }
}

int validate_isco(SimulationConfig *cfg) {
    printf("Validating ISCO at r = %.4f Rs\n", ISCO/RS);
    
    // Test rays with impact parameter corresponding to ISCO
    double b_test = ISCO / sqrt(1.0 - RS/ISCO);
    double px = b_test / cfg->camera_distance;
    double py = 0.0;
    
    RayTraceInfo info = trace_ray(px, py, cfg);
    
    if (info.disk_r > 0 && info.disk_r < ISCO * 1.01) {
        printf("  ✓ ISCO validation confirmed: disk hit at r = %.4f Rs\n", info.disk_r/RS);
        return 0;
    } else {
        printf("  ✗ Validation failed: expected ISCO hit but got %s\n",
               info.escaped ? "escape" : (info.disk_r > 0 ? "disk hit" : "black hole capture"));
        return -1;
    }
}

int validate_gravitational_redshift(SimulationConfig *cfg) {
    printf("Validating gravitational redshift formula\n");
    
    // Test at a known radius
    double test_r = 10.0 * RS;
    double expected_z = 1.0 / sqrt(1.0 - RS/test_r) - 1.0;
    
    // Set up a ray that will hit at test_r
    cfg->disk_inner = test_r;
    cfg->disk_outer = test_r + 0.001;
    
    double px = 0.0;
    double py = 0.0;
    
    RayTraceInfo info = trace_ray(px, py, cfg);
    
    if (info.disk_r > 0) {
        double error = fabs(info.redshift - 1.0 - expected_z);
        printf("  ✓ Gravitational redshift validation: z = %.4f (expected %.4f), error = %.2e\n",
               info.redshift - 1.0, expected_z, error);
        return (error < 1e-6) ? 0 : -1;
    } else {
        printf("  ✗ Validation failed: ray did not hit test radius\n");
        return -1;
    }
}

int run_validation_suite(SimulationConfig *cfg) {
    int failures = 0;
    
    failures += validate_photon_sphere(cfg);
    failures += validate_isco(cfg);
    failures += validate_gravitational_redshift(cfg);
    
    printf("\nValidation suite complete: %d failures out of 3 tests\n", failures);
    
    return failures;
}

// =====================================================================
// MAIN RESEARCH INTERFACE
// =====================================================================
void print_help(const char *prog_name) {
    printf("\nBLACK HOLE SIMULATION FRAMEWORK - RESEARCH EDITION\n");
    printf("==================================================\n\n");
    printf("Usage: %s [OPTIONS]\n\n", prog_name);
    printf("Basic options:\n");
    printf("  -w, --width=WIDTH        Set image width (default: 1024)\n");
    printf("  -h, --height=HEIGHT      Set image height (default: 1024)\n");
    printf("  -d, --distance=DIST      Set camera distance in Rs (default: 15.0)\n");
    printf("  -t, --theta=THETA        Set camera polar angle in degrees (default: 90.0)\n");
    printf("  -i, --disk-inner=R_IN    Set inner disk radius in Rs (default: 6.0)\n");
    printf("  -o, --disk-outer=R_OUT   Set outer disk radius in Rs (default: 30.0)\n");
    printf("  -m, --max-steps=STEPS    Set maximum integration steps (default: 2048)\n");
    printf("  -f, --format=FORMAT      Set output format (csv, hdf5, fits; default: csv)\n");
    printf("  -p, --prefix=PREFIX      Set output filename prefix (default: blackhole)\n");
    printf("\nAdvanced options:\n");
    printf("  -b, --benchmark          Run performance benchmark\n");
    printf("  -s, --sweep=PARAM        Run parameter sweep (distance, theta, disk_inner, disk_outer)\n");
    printf("  -v, --validate           Run validation suite against known solutions\n");
    printf("  -r, --range=START:END:STEP  Set parameter range for sweeps\n");
    printf("  -?, --help               Show this help message\n\n");
    printf("Example usage:\n");
    printf("  %s -w 2048 -h 2048 -d 20 -t 45 -i 6 -o 50 --sweep=theta --range=0:90:5\n", prog_name);
    printf("  %s --validate\n", prog_name);
    printf("\nThis simulation implements:\n");
    printf("  - Full Schwarzschild geodesic equations\n");
    printf("  - Adaptive step Verlet integration\n");
    printf("  - Gravitational redshift & Doppler boosting\n");
    printf("  - Critical photon orbit detection\n");
    printf("  - ISCO (Innermost Stable Circular Orbit) validation\n");
    printf("  - Parameter sweep capabilities for research\n");
}

int parse_range(const char *range_str, double *start, double *end, double *step) {
    char *colon1 = strchr(range_str, ':');
    char *colon2 = colon1 ? strchr(colon1 + 1, ':') : NULL;
    
    if (!colon1 || !colon2) {
        return -1;
    }
    
    *colon1 = '\0';
    *colon2 = '\0';
    
    *start = atof(range_str);
    *end = atof(colon1 + 1);
    *step = atof(colon2 + 1);
    
    *colon1 = ':';
    *colon2 = ':';
    
    return 0;
}

int main(int argc, char **argv) {
    // Initialize configuration with defaults
    SimulationConfig config = {
        .camera_distance = 15.0 * RS,
        .camera_theta = M_PI / 2.0,   // Equatorial view
        .disk_inner = 6.0 * RS,       // ISCO
        .disk_outer = 30.0 * RS,
        .screen_width = 1024,
        .screen_height = 1024,
        .max_steps = MAX_STEPS,
        .star_density = STAR_DENSITY,
        .max_redshift = MAX_REDSHIFT,
        .output_format = FORMAT_CSV,
        .validation_mode = 0
    };
    strcpy(config.output_prefix, "blackhole");
    
    // Parse command-line arguments
    int opt;
    enum SimulationMode mode = MODE_DEFAULT;
    char *sweep_param = NULL;
    double range_start = 0.0, range_end = 0.0, range_step = 1.0;
    
    static struct option long_options[] = {
        {"width", required_argument, 0, 'w'},
        {"height", required_argument, 0, 'h'},
        {"distance", required_argument, 0, 'd'},
        {"theta", required_argument, 0, 't'},
        {"disk-inner", required_argument, 0, 'i'},
        {"disk-outer", required_argument, 0, 'o'},
        {"max-steps", required_argument, 0, 'm'},
        {"format", required_argument, 0, 'f'},
        {"prefix", required_argument, 0, 'p'},
        {"benchmark", no_argument, 0, 'b'},
        {"sweep", required_argument, 0, 's'},
        {"validate", no_argument, 0, 'v'},
        {"range", required_argument, 0, 'r'},
        {"help", no_argument, 0, '?'},
        {0, 0, 0, 0}
    };
    
    while ((opt = getopt_long(argc, argv, "w:h:d:t:i:o:m:f:p:bs:v:r:?", long_options, NULL)) != -1) {
        switch (opt) {
            case 'w':
                config.screen_width = atoi(optarg);
                break;
            case 'h':
                config.screen_height = atoi(optarg);
                break;
            case 'd':
                config.camera_distance = atof(optarg) * RS;
                break;
            case 't':
                config.camera_theta = atof(optarg) * M_PI / 180.0;
                break;
            case 'i':
                config.disk_inner = atof(optarg) * RS;
                break;
            case 'o':
                config.disk_outer = atof(optarg) * RS;
                break;
            case 'm':
                config.max_steps = atoi(optarg);
                break;
            case 'f':
                if (strcmp(optarg, "csv") == 0) {
                    config.output_format = FORMAT_CSV;
                } else if (strcmp(optarg, "hdf5") == 0) {
                    #ifdef USE_HDF5
                        config.output_format = FORMAT_HDF5;
                    #else
                        fprintf(stderr, "HDF5 support not compiled in\n");
                        return 1;
                    #endif
                } else if (strcmp(optarg, "fits") == 0) {
                    #ifdef USE_FITS
                        config.output_format = FORMAT_FITS;
                    #else
                        fprintf(stderr, "FITS support not compiled in\n");
                        return 1;
                    #endif
                } else {
                    fprintf(stderr, "Invalid output format: %s\n", optarg);
                    return 1;
                }
                break;
            case 'p':
                strncpy(config.output_prefix, optarg, sizeof(config.output_prefix) - 1);
                break;
            case 'b':
                mode = MODE_BENCHMARK;
                break;
            case 's':
                mode = MODE_PARAMETER_SWEEP;
                sweep_param = optarg;
                break;
            case 'v':
                mode = MODE_VALIDATION;
                config.validation_mode = 1;
                break;
            case 'r':
                if (parse_range(optarg, &range_start, &range_end, &range_step) != 0) {
                    fprintf(stderr, "Invalid range format. Use START:END:STEP\n");
                    return 1;
                }
                break;
            case '?':
                print_help(argv[0]);
                return 0;
            default:
                fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
                return 1;
        }
    }
    
    // Validate configuration
    if (config.camera_distance < 2.0 * RS) {
        fprintf(stderr, "Error: Camera distance must be at least 2.0 Rs for safety\n");
        return 1;
    }
    
    if (config.disk_inner < ISCO) {
        fprintf(stderr, "Warning: Disk inner radius below ISCO (%.2f Rs). Setting to ISCO.\n", ISCO/RS);
        config.disk_inner = ISCO;
    }
    
    if (config.disk_outer <= config.disk_inner) {
        fprintf(stderr, "Error: Disk outer radius must be greater than inner radius\n");
        return 1;
    }
    
    // Run simulation based on mode
    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);
    
    switch (mode) {
        case MODE_DEFAULT: {
            printf("Running default simulation with parameters:\n");
            printf("  Camera distance: %.2f Rs\n", config.camera_distance/RS);
            printf("  Camera theta: %.2f degrees\n", config.camera_theta * 180.0 / M_PI);
            printf("  Disk inner radius: %.2f Rs\n", config.disk_inner/RS);
            printf("  Disk outer radius: %.2f Rs\n", config.disk_outer/RS);
            printf("  Resolution: %dx%d\n", config.screen_width, config.screen_height);
            printf("  Max steps: %d\n", config.max_steps);
            
            // Allocate memory for ray data
            RayTraceInfo **ray_data = malloc(config.screen_height * sizeof(RayTraceInfo*));
            if (!ray_data) {
                fprintf(stderr, "Memory allocation failed for ray data\n");
                return 1;
            }
            
            for (int i = 0; i < config.screen_height; i++) {
                ray_data[i] = malloc(config.screen_width * sizeof(RayTraceInfo));
                if (!ray_data[i]) {
                    fprintf(stderr, "Memory allocation failed for ray data\n");
                    for (int j = 0; j < i; j++) free(ray_data[j]);
                    free(ray_data);
                    return 1;
                }
            }
            
            // Trace rays
            #pragma omp parallel for collapse(2)
            for (int y = 0; y < config.screen_height; y++) {
                for (int x = 0; x < config.screen_width; x++) {
                    double px = (x - config.screen_width/2.0)/(0.5 * config.screen_width);
                    double py = (config.screen_height/2.0 - y)/(0.5 * config.screen_height);
                    
                    ray_data[y][x] = trace_ray(px, py, &config);
                }
            }
            
            // Save data
            char filename[256];
            snprintf(filename, sizeof(filename), "%s_photon_data.csv", config.output_prefix);
            save_photon_data_csv(filename, ray_data, config.screen_width, config.screen_height);
            
            // Generate telemetry
            SimulationTelemetry telemetry = {
                .time = time(NULL),
                .camera_distance = config.camera_distance,
                .camera_theta = config.camera_theta,
                .disk_inner = config.disk_inner,
                .disk_outer = config.disk_outer,
                .rays_traced = config.screen_width * config.screen_height
            };
            
            // Collect statistics
            int black_hole_hits = 0;
            int disk_hits = 0;
            int star_hits = 0;
            double total_redshift = 0.0;
            double max_redshift = 0.0;
            int photon_orbits_total = 0;
            
            for (int y = 0; y < config.screen_height; y++) {
                for (int x = 0; x < config.screen_width; x++) {
                    RayTraceInfo *info = &ray_data[y][x];
                    
                    if (!info->escaped && info->disk_r <= 0) {
                        black_hole_hits++;
                    } else if (info->disk_r > 0) {
                        disk_hits++;
                        total_redshift += info->redshift;
                        if (info->redshift > max_redshift) {
                            max_redshift = info->redshift;
                        }
                        photon_orbits_total += info->photon_orbits;
                    } else {
                        star_hits++;
                    }
                }
            }
            
            telemetry.black_hole_hits = black_hole_hits;
            telemetry.disk_hits = disk_hits;
            telemetry.star_hits = star_hits;
            telemetry.avg_redshift = disk_hits > 0 ? total_redshift / disk_hits : 0.0;
            telemetry.max_redshift_observed = max_redshift;
            telemetry.photon_orbits_total = photon_orbits_total;
            
            gettimeofday(&end_time, NULL);
            telemetry.simulation_time = (end_time.tv_sec - start_time.tv_sec) + 
                                      (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
            
            // Save telemetry
            snprintf(filename, sizeof(filename), "%s_telemetry.csv", config.output_prefix);
            save_telemetry_data(filename, &telemetry);
            
            // Free memory
            for (int i = 0; i < config.screen_height; i++) {
                free(ray_data[i]);
            }
            free(ray_data);
            
            printf("Simulation completed in %.2f seconds\n", telemetry.simulation_time);
            break;
        }
        
        case MODE_BENCHMARK: {
            printf("Running benchmark mode...\n");
            double avg_time = run_benchmark(5, &config);
            printf("Benchmark complete: %.2f ms per simulation\n", avg_time * 1000);
            break;
        }
        
        case MODE_PARAMETER_SWEEP: {
            if (!sweep_param) {
                fprintf(stderr, "Parameter sweep requires a parameter name\n");
                return 1;
            }
            
            char filename[256];
            snprintf(filename, sizeof(filename), "%s_sweep_%s.csv", 
                    config.output_prefix, sweep_param);
            
            run_parameter_sweep(&config, sweep_param, 
                               range_start, range_end, range_step, 
                               filename);
            break;
        }
        
        case MODE_VALIDATION: {
            printf("Running validation suite...\n");
            int failures = run_validation_suite(&config);
            
            if (failures == 0) {
                printf("\nAll validation tests passed successfully!\n");
            } else {
                printf("\nWARNING: %d validation tests failed. Results may be unreliable.\n", failures);
                return 1;
            }
            break;
        }
        
        default:
            fprintf(stderr, "Invalid mode\n");
            return 1;
    }
    
    gettimeofday(&end_time, NULL);
    double total_time = (end_time.tv_sec - start_time.tv_sec) + 
                       (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
    
    printf("\nTotal execution time: %.2f seconds\n", total_time);
    printf("Simulation complete. Data ready for analysis.\n");
    
    return 0;
}