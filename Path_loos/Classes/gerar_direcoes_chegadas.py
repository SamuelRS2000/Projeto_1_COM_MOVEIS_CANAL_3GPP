import os
import sys

# Garante a importacao a partir do diretorio real deste arquivo
_diretorio_real = os.path.dirname(os.path.abspath(__file__))
if _diretorio_real not in sys.path:
    sys.path.insert(0, _diretorio_real)

from gerar_atrasos import gerar_atrasos
import matplotlib.pyplot as plt
import numpy as np


class gerar_direcoes_chegadas(gerar_atrasos):

    def __init__(self, ambiente="UMi", freq=3, N_multi=100, **kwargs):
        # Heranca da classe de geracao de pequena escala
        super().__init__(
            ambiente=ambiente, freq=freq, N_multi=N_multi, **kwargs
        )

        # Angulos de visada de chegada convertidos para graus
        phi_LoS_prime_deg = np.degrees(self.phi_in)
        theta_LoS_prime_deg = np.degrees(self.theta_in)

        # --- Direcoes de Chegada - Angulo Azimutal (Slides 17 e 18) ---
        # Potencia normalizada travada entre 1e-15 e 1.0 para seguranca numerica no logaritmo
        rel_power = np.clip(self.alpha_sq / np.max(self.alpha_sq), 1e-15, 1.0)

        # Angulos azimutais iniciais (Eq. 4)
        phi_prime_prime = 1.42 * self.esp_Azt_in * np.sqrt(-np.log(rel_power))

        # Sinais discretos {-1, 1} e flutuacoes Gaussianas
        U_n_phi = np.random.choice([-1.0, 1.0], size=self.N_multi)
        Y_n_phi = np.random.normal(0.0, self.esp_Azt_in / 7.0, size=self.N_multi)

        # Angulos azimutais de chegada finais [graus] (Eq. 5)
        self.phi_n = U_n_phi * phi_prime_prime + Y_n_phi + phi_LoS_prime_deg

        # Conforme a TR 38.901 e a Eq. 8 do Slide 21 (theta'_1 = theta'_LoS), 
        # apenas a primeira componente (n=1) e forcada na direcao de visada direta.
        if self.flag_Loss == 1:
            self.phi_n[0] = phi_LoS_prime_deg

        # --- Direcoes de Chegada - Angulo em Elevacao (Slides 20 e 21) ---
        # Angulos de elevacao iniciais (Eq. 6)
        theta_prime_prime = -self.esp_elv_in * np.log(rel_power)

        # Sinais discretos {-1, 1} e flutuacoes Gaussianas
        U_n_theta = np.random.choice([-1.0, 1.0], size=self.N_multi)
        Y_n_theta = np.random.normal(
            0.0, self.esp_elv_in / 7.0, size=self.N_multi
        )

        # Angulos de elevacao de chegada finais [graus] (Eq. 7)
        self.theta_n = (
            U_n_theta * theta_prime_prime + Y_n_theta + theta_LoS_prime_deg
        )

        # Se houver visada direta, fixa a primeira componente no LoS (Eq. 8)
        if self.flag_Loss == 1:
            self.theta_n[0] = theta_LoS_prime_deg

        # --- Vetores de Direcao de Chegada r_n (Slide 16, Eq. 3) ---
        phi_rad = np.radians(self.phi_n)
        theta_rad = np.radians(self.theta_n)

        # Matriz 3 x N onde cada coluna e o vetor unitario r_n
        rx = np.cos(phi_rad) * np.sin(theta_rad)
        ry = np.sin(phi_rad) * np.sin(theta_rad)
        rz = np.cos(theta_rad)
        self.r_n = np.vstack([rx, ry, rz])

    # Espectro Angular de Potencia Azimutal (Slide 19 / Figura 2)
    def plot_Azt_chegadas(self, titulo=None, ax=None):
        show_plot = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": "polar"})
            show_plot = True

        phi_rad = np.radians(self.phi_n)

        # Tracado das componentes no plano polar
        for n in range(self.N_multi):
            ax.plot(
                [phi_rad[n], phi_rad[n]],
                [1e-5, self.alpha_sq[n]],
                color="black",
                linewidth=1.0,
            )
            ax.scatter(phi_rad[n], self.alpha_sq[n], color="black", s=20)

        # Escala logaritmica radial
        ax.set_rscale("log")
        ax.set_rlim(1e-5, max(0.1, np.max(self.alpha_sq) * 1.5))

        ax.set_theta_zero_location("E")
        
        titulo_padrao = titulo if titulo else "Espectro Angular de Potência (Azimute)\n"
        ax.set_title(titulo_padrao, fontsize=12)

        # Caixa com espalhamento angular azimutal
        ax.text(
            0.05,
            0.05,
            rf"$\sigma_{{\phi,AoA}} = {self.esp_Azt_in:.2f}^\circ$",
            transform=ax.transAxes,
            fontsize=11,
            bbox=dict(
                boxstyle="square,pad=0.3", facecolor="white", edgecolor="black"
            ),
        )

        if show_plot:
            plt.tight_layout()
            plt.show()
            
        return ax

    # Espectro Angular de Potencia em Elevacao (Slide 22 / Figura 3)
    def plot_elv_chegadas(self, titulo=None, ax=None):
        show_plot = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6), subplot_kw={"projection": "polar"})
            show_plot = True

        theta_rad = np.radians(self.theta_n)

        # Tracado das componentes no plano polar
        for n in range(self.N_multi):
            ax.plot(
                [theta_rad[n], theta_rad[n]],
                [1e-5, self.alpha_sq[n]],
                color="black",
                linewidth=1.0,
            )
            ax.scatter(theta_rad[n], self.alpha_sq[n], color="black", s=20)

        # Escala logaritmica radial
        ax.set_rscale("log")
        ax.set_rlim(1e-5, max(0.1, np.max(self.alpha_sq) * 1.5))

        ax.set_theta_zero_location("E")
        
        titulo_padrao = titulo if titulo else "Espectro Angular de Potência (Elevação)\n"
        ax.set_title(titulo_padrao, fontsize=12)

        # Caixa com espalhamento angular de elevacao
        ax.text(
            0.05,
            0.05,
            rf"$\sigma_{{\theta,AoA}} = {self.esp_elv_in:.2f}^\circ$",
            transform=ax.transAxes,
            fontsize=11,
            bbox=dict(
                boxstyle="square,pad=0.3", facecolor="white", edgecolor="black"
            ),
        )

        if show_plot:
            plt.tight_layout()
            plt.show()
            
        return ax

    # Vetores Direcao das Componentes Multipercurso r_n (Slide 24 / Figura 4)
    def plot_r_n(self, titulo=None, ax=None):
        show_plot = False
        if ax is None:
            fig = plt.figure(figsize=(7, 7))
            ax = fig.add_subplot(projection="3d")
            show_plot = True

        # Tracado das componentes dispersas (vermelho com marcadores triangulares)
        start_idx = 1 if self.flag_Loss == 1 else 0
        for n in range(start_idx, self.N_multi):
            vx, vy, vz = self.r_n[0, n], self.r_n[1, n], self.r_n[2, n]
            ax.plot([0, vx], [0, vy], [0, vz], color="red", linewidth=1.0)
            ax.scatter(
                vx,
                vy,
                vz,
                color="red",
                marker="^",
                s=40,
                facecolors="none",
                edgecolors="red",
            )

        # Componente principal LoS destacada em azul caso exista
        if self.flag_Loss == 1:
            vx, vy, vz = self.r_n[0, 0], self.r_n[1, 0], self.r_n[2, 0]
            ax.plot([0, vx], [0, vy], [0, vz], color="blue", linewidth=1.5)
            ax.scatter(
                vx,
                vy,
                vz,
                color="blue",
                marker="^",
                s=50,
                facecolors="none",
                edgecolors="blue",
            )

        ax.set_xlim([-1, 1])
        ax.set_ylim([-1, 1])
        ax.set_zlim([-1, 1])

        ax.set_xlabel("Eixo X")
        ax.set_ylabel("Eixo Y")
        ax.set_zlabel("Eixo Z")
        
        titulo_padrao = titulo if titulo else "Vetores Direção das Componentes Multipercurso\n"
        ax.set_title(titulo_padrao, fontsize=12)

        if show_plot:
            plt.tight_layout()
            plt.show()
            
        return ax


# --- Teste isolado do arquivo gerar_direcoes_chegadas.py ---
if __name__ == "__main__":
    canal = gerar_direcoes_chegadas(ambiente="UMi", freq=3.5, N_multi=100)

    print("\n--- Teste Direções de Chegada (UMi) ---")
    print(f"Condição: {'LoS' if canal.flag_Loss == 1 else 'NLoS'}")
    print(
        f"Espalhamento Azimute Chegada: {canal.esp_Azt_in:.2f}° | Elevação: {canal.esp_elv_in:.2f}°"
    )
    print(
        f"Primeiro Azimute (phi'_1): {canal.phi_n[0]:.2f}° | Elevação (theta'_1): {canal.theta_n[0]:.2f}°"
    )
    print(
        f"Norma média dos vetores r_n: {np.mean(np.linalg.norm(canal.r_n, axis=0)):.4f}"
    )

    # Criacao do painel unificado para os tres graficos
    fig = plt.figure(figsize=(12, 10))
    
    # Eixos superiores polares (linha 1, colunas 1 e 2)
    ax1 = plt.subplot(2, 2, 1, projection="polar")
    ax2 = plt.subplot(2, 2, 2, projection="polar")
    
    # Eixo inferior 3D (ocupa posicoes 3 e 4 da malha)
    ax3 = plt.subplot(2, 2, (3, 4), projection="3d")
    
    # Injecao dos graficos em cada sub-eixo
    canal.plot_Azt_chegadas(titulo="Espectro Angular (Azimute)\n", ax=ax1)
    canal.plot_elv_chegadas(titulo="Espectro Angular (Elevação)\n", ax=ax2)
    canal.plot_r_n(titulo="Vetores Direção Tridimensionais (r_n)\n", ax=ax3)

    plt.tight_layout()
    plt.show()