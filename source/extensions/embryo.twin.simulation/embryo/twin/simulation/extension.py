import omni.ext
import omni.ui as ui
import omni.usd
import omni.kit.app

from pxr import UsdGeom
from pxr import Gf

from .embryo import EmbryoSimulation


class EmbryoTwinExtension(omni.ext.IExt):

    def on_startup(self, ext_id):

        print(
            "[EmbryoTwin] Starting Embryo Digital Twin"
        )

        # ----------------------------------------
        # Simulation
        # ----------------------------------------

        self.simulation = EmbryoSimulation()

        self.running = False

        # Simulation speed:
        # biological hours per real second.
        self.simulation_speed = 1.0

        # ----------------------------------------
        # USD
        # ----------------------------------------

        self.usd_context = omni.usd.get_context()

        self.stage = None

        # ----------------------------------------
        # UI
        # ----------------------------------------

        self.window = ui.Window(
            "Embryo Digital Twin",
            width=360,
            height=420
        )

        self._build_ui()

        # ----------------------------------------
        # Kit update callback
        # ----------------------------------------

        update_stream = (
            omni.kit.app
            .get_app()
            .get_update_event_stream()
        )

        self._update_subscription = (
            update_stream.create_subscription_to_pop(
                self._on_update,
                name="Embryo Twin Update"
            )
        )

        # Create initial scene

        self._create_stage()

        self._sync_usd()

    # =====================================================
    # UI
    # =====================================================

    def _build_ui(self):

        with self.window.frame:

            with ui.VStack(
                spacing=8,
                height=0
            ):

                ui.Label(
                    "EMBRYO DIGITAL TWIN",
                    height=30,
                    style={
                        "font_size": 20
                    }
                )

                ui.Separator()

                self.time_label = ui.Label(
                    "Simulation time: 0.00 h"
                )

                self.cell_label = ui.Label(
                    "Cells: 1"
                )

                self.status_label = ui.Label(
                    "Status: STOPPED"
                )

                ui.Separator()

                with ui.HStack(
                    spacing=5,
                    height=30
                ):

                    ui.Button(
                        "START",
                        clicked_fn=self._start
                    )

                    ui.Button(
                        "STOP",
                        clicked_fn=self._stop
                    )

                    ui.Button(
                        "RESET",
                        clicked_fn=self._reset
                    )

                ui.Separator()

                ui.Label(
                    "Simulation speed (hours / real second)"
                )

                self.speed_model = ui.SimpleFloatModel(
                    self.simulation_speed
                )

                self.speed_slider = ui.FloatSlider(
                    model=self.speed_model,
                    min=0.1,
                    max=10.0
                )

                self.speed_value = ui.Label(
                    f"{self.simulation_speed:.1f}"
                )

                self.speed_model.add_value_changed_fn(
                    self._speed_changed
                )

                ui.Separator()

                ui.Label(
                    "Cell parameters"
                )

                self.growth_label = ui.Label(
                    f"Growth rate: "
                    f"{self.simulation.growth_rate:.3f}"
                )

                self.division_label = ui.Label(
                    f"Division radius: "
                    f"{self.simulation.division_radius:.3f}"
                )

                ui.Separator()

                ui.Label(
                    "Biology: growth → mitosis → daughter cells"
                )

    # =====================================================
    # UI CALLBACKS
    # =====================================================

    def _start(self):

        self.running = True

        self.status_label.text = (
            "Status: RUNNING"
        )

        print("[EmbryoTwin] Simulation started")

    def _stop(self):

        self.running = False

        self.status_label.text = (
            "Status: STOPPED"
        )

        print("[EmbryoTwin] Simulation stopped")

    def _reset(self):

        self.running = False

        self.simulation.reset()

        self._create_stage()

        self._sync_usd()

        self.status_label.text = (
            "Status: RESET"
        )

        print("[EmbryoTwin] Simulation reset")

    def _speed_changed(self, model):

        self.simulation_speed = (
            model.as_float
        )

        self.speed_value.text = (
            f"{self.simulation_speed:.1f}"
        )

    # =====================================================
    # UPDATE LOOP
    # =====================================================

    def _on_update(self, event):

        if not self.running:

            return

        dt_real = event.payload["dt"]

        # Avoid a huge simulation jump if
        # the application stalls.
        dt_real = min(
            float(dt_real),
            0.1
        )

        biological_dt = (
            dt_real *
            self.simulation_speed
        )

        self.simulation.update(
            biological_dt
        )

        self._sync_usd()

        self.time_label.text = (
            f"Simulation time: "
            f"{self.simulation.time:.2f} h"
        )

        self.cell_label.text = (
            f"Cells: "
            f"{len(self.simulation.alive_cells())}"
        )

    # =====================================================
    # USD
    # =====================================================

    def _create_stage(self):

        self.usd_context.new_stage()

        self.stage = (
            self.usd_context.get_stage()
        )

        # Root

        UsdGeom.Xform.Define(
            self.stage,
            "/World"
        )

        UsdGeom.Xform.Define(
            self.stage,
            "/World/Embryo"
        )

        # Camera-friendly ground/reference

        self.stage.SetDefaultPrim(
            self.stage.GetPrimAtPath(
                "/World"
            )
        )

    def _sync_usd(self):

        if self.stage is None:

            return

        embryo_path = "/World/Embryo"

        # Existing cell paths

        existing_paths = set()

        embryo_prim = self.stage.GetPrimAtPath(
            embryo_path
        )

        if embryo_prim:

            for child in embryo_prim.GetChildren():

                existing_paths.add(
                    str(child.GetPath())
                )

        active_paths = set()

        # ----------------------------------------
        # Create/update cells
        # ----------------------------------------

        for cell in self.simulation.alive_cells():

            path = (
                f"{embryo_path}/"
                f"Cell_{cell.cell_id:06d}"
            )

            active_paths.add(path)

            sphere = UsdGeom.Sphere.Define(
                self.stage,
                path
            )

            sphere.GetRadiusAttr().Set(
                float(cell.radius)
            )

            # --------------------------------
            # Cell color based on fate
            # --------------------------------

            if cell.fate == "PROGENITOR":

                color = Gf.Vec3f(
                    0.1, 0.3, 1.0
                )

            elif cell.fate == "TYPE_A":

                color = Gf.Vec3f(
                    1.0, 0.1, 0.1
                )

            elif cell.fate == "TYPE_B":

                color = Gf.Vec3f(
                    0.1, 1.0, 0.2
                )

            else:

                color = Gf.Vec3f(
                    0.7, 0.7, 0.7
                )

            """

            # TEST COLORS
            if cell.generation % 3 == 0:

                color = Gf.Vec3f(0.1, 0.3, 1.0)   # blue

            elif cell.generation % 3 == 1:

                color = Gf.Vec3f(1.0, 0.1, 0.1)   # red

            else:

                color = Gf.Vec3f(0.1, 1.0, 0.2)   # green

            """

            sphere.GetDisplayColorAttr().Set(
                [color]
            )


            # --------------------------------
            # Position
            # --------------------------------

            xformable = UsdGeom.Xformable(
                sphere.GetPrim()
            )

            translate_op = None

            for op in xformable.GetOrderedXformOps():

                if op.GetOpType() == UsdGeom.XformOp.TypeTranslate:

                    translate_op = op
                    break

            if translate_op is None:

                translate_op = xformable.AddTranslateOp()

            translate_op.Set(
                Gf.Vec3d(
                    float(cell.position[0]),
                    float(cell.position[1]),
                    float(cell.position[2])
                )
            )

            prim = sphere.GetPrim()

            # --------------------------------
            # Biological metadata
            # --------------------------------

            custom_data = {
                "cell_id": int(cell.cell_id),

                "parent_id": (
                    int(cell.parent_id)
                    if cell.parent_id is not None
                    else -1
                ),

                "generation": int(cell.generation),

                "age_hours": float(cell.age),

                "oxygen": float(cell.oxygen),

                "morphogen_A": float(cell.morphogen_a),

                "cell_cycle": str(cell.cell_cycle),

                "fate": str(cell.fate),

                "lineage_id": str(cell.lineage_id),
            }

            prim.SetCustomData(custom_data)
        # ----------------------------------------
        # Remove cells that divided
        # ----------------------------------------

        for path in existing_paths - active_paths:

            prim = self.stage.GetPrimAtPath(
                path
            )

            if prim.IsValid():

                self.stage.RemovePrim(
                    path
                )

    # =====================================================
    # SHUTDOWN
    # =====================================================

    def on_shutdown(self):

        print(
            "[EmbryoTwin] Shutting down"
        )

        self._update_subscription = None

        if self.window:

            self.window.destroy()

            self.window = None