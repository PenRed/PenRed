## Generic Simulation Parameters

Simulation parameters are specific to the main configuration of the simulation, unlike source, tally, geometry, or material configurations. This section describes the available parameters for the `pen_main` program, which control the overall behavior and execution of the simulation. This parameters can be found in the `Simulation properties` world tag.

<img src="/simulation-configuration/images/worldPropertiesSimParameters.png" alt="Simulation properties tag" width="400" style="display: block; margin: 0 auto"/>

### Dumps and Results Control

PenRed provides the capability to dump the current state of the entire simulation. This dump can be used to:

- Resume a simulation after computer crash.
- Store the final results in binary format.
- Combine the results of simulations with the same configuration.

To configure this feature, the following parameters are available:

#### Dumps

This parameters group handles the simulation checkpointing. By default, both, writing and reading dump files is disabled.

- **`Write interval`**  
  Specifies the time interval, in seconds, between successive dumps.

- **`Write File`**  
  Specifies the filename used to save the dumps. Each thread will automatically append the `thN` prefix to the filename, where `N` is the thread number.
  
- **`Read File`**  
  Expects a string with the name of a dump file to resume a previous simulation.  
  **Note:** If multiple threads are used, the simulation must be resumed with the same number of threads as the original simulation. Each thread will attempt to read a file with the specified name prefixed by `thN`, where `N` is the thread number. For example, if the filename is `dump.dat`, the file for thread number 2 must be named `th2dump.dat`.
  
#### Results

The *Results* group controls how the simulation results are saved.

- **`Prefix`**  
  This parameter allows the user to specify a prefix for the ASCII results produced by every tally. The prefix will be appended to the usual results filenames.  
  **Note:** To specify an output directory, a trailing slash is required:
    + Use `/` for Linux (e.g., `output/`).
    + Use `\` for Windows (e.g., `output\`).  
  The specified directory must exist before the simulation starts.

- **`Final Dump`**  
  If enabled, forces a dump in each thread at the end of the simulation.  
  **Default:** Disabled.

- **`ASCII Results`**  
  Enables or disables result reports in ASCII format. If disabled, the result reports will be written to a dump file named `results.dump`. Additionally, no partial results will be printed during the simulation.  
  **Default:** Enabled.

- **`ASCII Partial Results`**  
  If `ASCII Results` is enabled, this boolean value can enable or disable the ASCII printing of partial results.  
  **Default:** Disabled.

### Threading

The following parameters control the multi-threading capabilities of the simulation:

- **`Threads Selection`**  
  Controls how the number of threads are selected:
    + **`Auto`**: The number of threads is authomatically set according to the number of virtual cores.
    + **`Manual`**: The number of threads is specified by the user.
- **`Threads`**  
  An integer value that specifies the number of threads to be used. It is only available if the thread selection is set to *Manual*.  
- **`Seed Pair`**  
  Specify the seed pair number used to initialize the random number generator.

The number of histories to simulate on each source will be distributed among all specified threads.

### Time Limit

If enabled the user can specify the maximum duration, in seconds, that the simulation is allowed to run.

