use ui::MyEguiApp;

mod constants;
mod math;
mod ui;

fn main() {
    let native_options = eframe::NativeOptions::default();
    eframe::run_native(
        "My egui App",
        native_options,
        Box::new(|cc| Box::new(MyEguiApp::new(cc, 120.0, 100.0))),
    );
}
