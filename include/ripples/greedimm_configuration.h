
namespace ripples {
struct GreedimmConfiguration : public IMMConfiguration {
  void addCmdOptions(CLI::App &app) {
    IMMConfiguration::addCmdOptions(app);
    app.add_option("--seed-select-max-gpu-workers", seed_select_max_gpu_workers,
                   "The max number of GPU workers for seed selection.")
        ->group("Streaming-Engine Options");
  }
}
}