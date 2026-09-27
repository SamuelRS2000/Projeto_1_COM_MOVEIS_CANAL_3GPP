import matplotlib.pyplot as plt
import numpy as np


class gerar_param:
    """
    Entradas :  ambiente -> Tipo do ambiente
                freq -> Frequência central [GHz]
                x_UE, y_UE -> Coordenadas cartesianas do UE [m] (opcionais)
                vet_v_tx -> Vetor unitário de direção da velocidade do UE (opcional)

    Saídas :
    Parâmetros calculados disponíveis para as próximas etapas:
                h_BS -> Altura da BS [m]
                h_UT -> Altura do UT [m]
                d_BS_UE -> Distância horizontal 2D entre BS e UE [m]
                d_3D -> Distância tridimensional total entre BS e UE [m]
                mod_v_tx -> Módulo da velocidade do Móvel [km/h]
                vet_v_tx -> Vetor Unitário da Velocidade do Móvel [x,y,z]
                flag_Loss -> Se possui ou não Linha de visada (1 para LOS, 0 para NLOS)
                phi_out -> Azimute de Saída
                phi_in -> Azimute de Chegada
                theta_out -> Elevação de Saída
                theta_in -> Elevação de Chegada
    """

    def __init__(
        self,
        ambiente="UMi",
        freq=3,
        mod_v_tx=3,
        x_UE=None,
        y_UE=None,
        vet_v_tx=None,
    ):
        self.ambiente = ambiente
        self.freq = freq  # [GHz]
        self.mod_v_tx = (
            mod_v_tx  # [km/h] - fixado em 3 km/h 
        )

        # Definição do vetor unitário de velocidade do UT (Considerando arbitrarimente que o UT esta se movendo apenas no eixo X)
        if vet_v_tx is None:
            self.vet_v_tx = np.array([1.0, 0.0, 0.0])
        else:
            vet_v_tx = np.array(vet_v_tx, dtype=float)
            self.vet_v_tx = vet_v_tx / np.linalg.norm(vet_v_tx)

        # Conforme referenciado na TR 38.901, a TR 36.873 define: h_UT=3(nfl - 1) + 1.5
        # h_UT = 1.5 m.
        self.h_UT = 1.5  # [m]

        # Obtenção do Pr_LOS para cada ambiente conforme Tabelas 7.2-2 e 7.4.2-1
        # --- UMi - Street canyon
        if self.ambiente == "UMi":
            self.h_BS = 10  # [m]
            self.min_d_BS_UE = 10  # [m]
            self.ISD = 200

            # Definição das coordenadas com relação à origem na posição da BS
            self.posicionar_UE(
                self.min_d_BS_UE, self.ISD, x_in=x_UE, y_in=y_UE
            )

            self.Pr_LOS = np.where(
                self.d_BS_UE <= 18,
                1.0,
                (18 / self.d_BS_UE)
                + np.exp(-self.d_BS_UE / 36) * (1 - 18 / self.d_BS_UE),
            )

        # --- UMa
        elif self.ambiente == "UMa":
            self.h_BS = 25
            self.min_d_BS_UE = 35
            self.ISD = 500

            # Definição das coordenadas com relação à origem na posição da BS
            self.posicionar_UE(
                self.min_d_BS_UE, self.ISD, x_in=x_UE, y_in=y_UE
            )

            c_prime = np.where(
                self.h_UT <= 13, 0.0, ((self.h_UT - 13) / 10) ** 1.5
            )

            self.Pr_LOS = np.where(
                self.d_BS_UE <= 18,
                1.0,
                (
                    (18 / self.d_BS_UE)
                    + np.exp(-self.d_BS_UE / 63) * (1 - 18 / self.d_BS_UE)
                )
                * (
                    1
                    + c_prime
                    * (5 / 4)
                    * ((self.d_BS_UE / 100) ** 3)
                    * np.exp(-self.d_BS_UE / 150)
                ),
            )

        # --- Indoor - Mixed office
        elif self.ambiente == "Indoor - Mixed office":
            self.h_UT = 1
            self.h_BS = 3
            self.min_d_BS_UE = 0
            self.ISD = 20

            # Definição das coordenadas com relação à origem na posição da BS
            self.posicionar_UE(
                self.min_d_BS_UE, self.ISD, x_in=x_UE, y_in=y_UE
            )

            self.Pr_LOS = np.where(
                self.d_BS_UE <= 1.2,
                1.0,
                np.where(
                    self.d_BS_UE < 6.5,
                    np.exp(-(self.d_BS_UE - 1.2) / 4.7),
                    np.exp(-(self.d_BS_UE - 6.5) / 32.6) * 0.32,
                ),
            )

        # Conversão de array para float escalar, caso aplicável
        if isinstance(self.Pr_LOS, np.ndarray):
            self.Pr_LOS = float(self.Pr_LOS.item())

        # Sorteio de Bernoulli para determinar a existência de visada direta
        self.flag_Loss = int(np.random.binomial(1, self.Pr_LOS))

    def posicionar_UE(self, min_d, max_d, x_in=None, y_in=None):
        # Determinação das coordenadas X e Y do UE
        if x_in is not None and y_in is not None:
            self.x = float(x_in)
            self.y = float(y_in)
            self.d_BS_UE = np.sqrt(self.x**2 + self.y**2)
        else:
            # Geração aleatória uniforme dentro dos limites do site
            self.d_BS_UE = float(np.random.uniform(min_d, max_d))
            phi_sorteio = np.random.uniform(0, 2 * np.pi)
            self.x = self.d_BS_UE * np.cos(phi_sorteio)
            self.y = self.d_BS_UE * np.sin(phi_sorteio)

        # Cota Z relativa do UE em relação à BS
        self.z = float(self.h_UT - self.h_BS)

        # Distância tridimensional (3D) real
        self.d_3D = float(np.sqrt(self.d_BS_UE**2 + self.z**2))

        # Determinação dos ângulos de visada direta
        self.phi_out = float(np.arctan2(self.y, self.x))
        self.phi_in = float(np.arctan2(-self.y, -self.x))
        self.theta_out = float(np.arctan2(self.z, self.d_BS_UE))
        self.theta_in = float(-self.theta_out)

    def draw_UE_3d(self, ax=None):
        """
        Gera um gráfico 3D representando o enlace, mantendo o sistema de coordenadas 
        referenciado com a Estação Rádio Base (BS) na origem (0,0,0).
        """
        if ax is None:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(projection="3d")
            
            # Configuração dos eixos e título
            ax.set_title(f"Geometria do Enlace 3D ({self.ambiente})\n", fontsize=14, fontweight='bold')
            ax.set_xlabel("Eixo X [m]", labelpad=10)
            ax.set_ylabel("Eixo Y [m]", labelpad=10)
            ax.set_zlabel("Eixo Z [m]", labelpad=10)
            
            # Plotagem da BS fixa na origem (0,0,0) apenas uma vez
            ax.scatter(0, 0, 0, color="black", s=60, marker="o", zorder=5)
            ax.text(0, 0, 1.5, " BS (0,0,0)", color="black", fontweight="bold")

        # Definição das cores baseada na probabilidade LoS/NLoS
        if self.flag_Loss == 1:
            cor_enlace = "#2ca02c"  # Verde
            estilo_enlace = "-"
            lbl_enlace = "Visada Direta (LoS)"
        else:
            cor_enlace = "#ff7f0e"  # Laranja
            estilo_enlace = "--"
            lbl_enlace = "Obstruído (NLoS)"

        # Conexão entre a BS (origem) e as coordenadas relativas do UE
        ax.plot(
            [0, self.x], [0, self.y], [0, self.z], 
            color=cor_enlace, linestyle=estilo_enlace, linewidth=1.8, label=lbl_enlace
        )

        # Plotagem do marcador do UE
        ax.scatter(self.x, self.y, self.z, color="#1f77b4", s=50, marker="o", zorder=5)
        ax.text(self.x, self.y, self.z - 1.5, " UE", color="black", fontweight="bold")

        # Rótulo de distância no ponto médio
        ax.text(
            self.x / 2, self.y / 2, self.z / 2, 
            f" d={self.d_3D:.1f}m", color="black", fontweight="bold", ha='center',
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7)
        )

        # Organiza as legendas para não duplicar em testes com múltiplos UEs
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys(), loc="upper right", framealpha=0.9)

        return ax


# --- Teste isolado do arquivo gerar_param.py ---
if __name__ == "__main__":
    
    # Canal 1: Coordenadas especificadas manualmente
    canal1 = gerar_param(ambiente="UMi", freq=3.5, x_UE=25.0, y_UE=30.0)
    print("\n--- Teste gerar_param (Canal 1 - Coordenadas Manuais) ---")
    print(f"Posição UE: x={canal1.x:.1f} m, y={canal1.y:.1f} m, z={canal1.z:.1f} m")
    print(f"Distância 2D: {canal1.d_BS_UE:.2f} m | Distância 3D: {canal1.d_3D:.2f} m")
    print(f"Probabilidade LOS: {canal1.Pr_LOS:.4f} | LOS flag: {canal1.flag_Loss}")
    ax = canal1.draw_UE_3d()

    # Canal 2: Sorteio automático de posicionamento
    canal2 = gerar_param(ambiente="UMi", freq=3.5)
    print("\n--- Teste gerar_param (Canal 2 - Sorteio Interno) ---")
    print(f"Posição UE: x={canal2.x:.1f} m, y={canal2.y:.1f} m, z={canal2.z:.1f} m")
    print(f"Distância 2D: {canal2.d_BS_UE:.2f} m | Distância 3D: {canal2.d_3D:.2f} m")
    print(f"Probabilidade LOS: {canal2.Pr_LOS:.4f} | LOS flag: {canal2.flag_Loss}")
    canal2.draw_UE_3d(ax=ax)

    # Canal 3: Sorteio automático de posicionamento
    canal3 = gerar_param(ambiente="UMi", freq=3.5)
    print("\n--- Teste gerar_param (Canal 3 - Sorteio Interno) ---")
    print(f"Posição UE: x={canal3.x:.1f} m, y={canal3.y:.1f} m, z={canal3.z:.1f} m")
    print(f"Distância 2D: {canal3.d_BS_UE:.2f} m | Distância 3D: {canal3.d_3D:.2f} m")
    print(f"Probabilidade LOS: {canal3.Pr_LOS:.4f} | LOS flag: {canal3.flag_Loss}")
    canal3.draw_UE_3d(ax=ax)

    plt.show()