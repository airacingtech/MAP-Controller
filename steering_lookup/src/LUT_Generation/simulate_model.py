from dynamics.vehicle_dynamics_stown import vehicle_dynamics_st
from helpers.load_model import get_dotdict
from scipy.integrate import odeint
import numpy as np

model_name = "AV24_pacejka"
model, tiretype = model_name.split("_")

model = get_dotdict(model_name)

class Simulator:
    def func_ST(self, x, t, u, p):
        f = vehicle_dynamics_st(x, u, p, tiretype)
        return f

    def generate_lookup(self):
        start_steer = 0.0
        steer_fine_end = 0.1
        end_steer = 0.1
        n_steps_steer = 30  # used for fine region and then again for coarse region
        start_vel = 0.5
        end_vel = 70.0
        n_steps_vel = 65
        start_bank = 0.0  # Range for bank angle in radians (slight negative tilt)
        end_bank = -0.35  # Slight positive tilt
        n_steps_bank = 15  # Number of steps for bank angle
        # end_bank = 0.0
        # n_steps_bank = 1
        start_curv = 0.0  # Range for bank angle in radians (slight negative tilt)
        # end_curv = 0.0
        # n_steps_curv = 1
        end_curv = 0.006  # Slight positive tilt
        n_steps_curv = 10  # Number of steps for curv

        fine_steers = np.linspace(start_steer, steer_fine_end, n_steps_steer)
        steers = np.linspace(steer_fine_end, end_steer, n_steps_steer)
        steers = np.concatenate((fine_steers, steers))
        vels = np.linspace(start_vel, end_vel, n_steps_vel)
        banks = np.linspace(start_bank, end_bank, n_steps_bank)
        curvatures = np.linspace(start_curv, end_curv, n_steps_curv)

        n_steps_steer = len(steers)
        n_steps_bank = len(banks)
        self.lookup_table = np.empty([n_steps_steer + 1, n_steps_vel + 1, n_steps_bank + 1, n_steps_curv +1])

        self.lookup_table[0, 1:, 0,0] = vels  # Velocity row
        self.lookup_table[1:, 0, 0,0] = steers  # Steering column
        self.lookup_table[0, 0, 1:,0] = banks  # Bank angle row
        self.lookup_table[0, 0, 0,1:] = curvatures  # Bank angle row


        for steer_idx, steer in enumerate(steers):
            for vel_idx, vel in enumerate(vels):
                for bank_idx, bank in enumerate(banks):
                    for k_idx, k in enumerate(curvatures):
                      initialState = [0, 0, 0, vel, 0, 0]
                      u = [steer, 0, bank, k]  # Now includes bank as part of input
                      
                      dt = 0.01
                      duration = 2.0
                      t = np.arange(0, duration, dt)

                      self.sol = odeint(self.func_ST, initialState, t, args=(u, model))
                      # print(steer)
                      # print(vel)
                      # print(bank)
                      # print(k)
                      if abs(self.sol[-1, 5] - self.sol[-10, 5]) > 0.05:
                          # print("failed to converge, limit")
                          self.lookup_table[steer_idx+1, vel_idx+1, bank_idx+1,k_idx +1] = None
                          break

                      else:
                          a_lat = self.sol[-1, 5] * vel
                          a_lat_all = self.sol[:, 5] * vel
                          # print("saved")

                          self.lookup_table[steer_idx+1, vel_idx+1, bank_idx+1, k_idx+1] = a_lat

    def save_lookup(self):
      model, tires = model_name.split("_")
      file_path = "./models/" + model + "/" + model_name + "_lookup_table.csv"
      
      # Flatten the arrays for steer, vel, bank, curv, and lateral_acc
      steer_vals, vel_vals, bank_vals, curv_vals, lateral_acc_vals = [], [], [], [], []
      
      for steer_idx, steer in enumerate(self.lookup_table[1:, 0, 0, 0]):
          for vel_idx, vel in enumerate(self.lookup_table[0, 1:, 0, 0]):
              for bank_idx, bank in enumerate(self.lookup_table[0, 0, 1:, 0]):
                  for curv_idx, curv in enumerate(self.lookup_table[0, 0, 0, 1:]):
                      lateral_acc = self.lookup_table[steer_idx+1, vel_idx+1, bank_idx+1, curv_idx+1]
                      steer_vals.append(steer)
                      vel_vals.append(vel)
                      bank_vals.append(bank)
                      curv_vals.append(curv)
                      lateral_acc_vals.append(lateral_acc)

      # Stack the columns together
      data = np.column_stack((steer_vals, vel_vals, bank_vals, curv_vals, lateral_acc_vals))
      
      # Save to CSV
      np.savetxt(file_path, data, delimiter=",", header="steer,vel,bank,curv,lateral_acc", comments="")

if __name__ == "__main__":
    sim = Simulator()
    sim.generate_lookup()
    sim.save_lookup()
    print("Done")
