import os
import sys

# Diretório base dinâmico
_diretorio_raiz = os.path.dirname(os.path.abspath(__file__))

# Caminhos relativos 
_pasta_path_loss = os.path.join(_diretorio_raiz, "Path_loos")
_pasta_classes = os.path.join(_diretorio_raiz, "Classes")

# Adiciona ao sys.path de forma portátil
for pasta in [_diretorio_raiz, _pasta_path_loss, _pasta_classes]:
    if os.path.exists(pasta) and pasta not in sys.path:
        sys.path.insert(0, pasta)

from path_loss_3GPP_TR_38_901 import Canal3GPP
import matplotlib.pyplot as plt
import numpy as np


class P_RX(Canal3GPP):

    def __init__(self, ambiente="UMi", freq=3.5, N_multi=100, override_espa_atraso=None, **kwargs):
        # O kwarg override_espa_atraso é repassado pela herança até chegar em gerar_atrasos
        super().__init__(
            ambiente=ambiente, freq=freq, N_multi=N_multi, 
            override_espa_atraso=override_espa_atraso, **kwargs
        )

    # Gera o pulso retangular transmitido s_tilde(t) (Slides 28, 29 e 31)
    def gerar_pulso_tx(self, delta_t, atraso_ini=0.0, Nt=50000):
        t = np.linspace(0.0, 5.0 * delta_t, Nt)
        s_t = np.where(
            (t >= atraso_ini) & (t <= (atraso_ini + delta_t)), 1.0, 0.0
        )
        return t, s_t

    # Calcula o sinal recebido r_tilde(t) (Slide 30 - Eq. 12)
    def simular_sinal_rx(self, delta_t, atraso_ini=0.0, Nt=50000):
        t, s_t = self.gerar_pulso_tx(delta_t, atraso_ini=atraso_ini, Nt=Nt)

        # Amplitude da envoltória de cada componente multipercurso
        alpha_n = np.sqrt(self.alpha_sq)

        # Acumulador do sinal recebido complexo em banda básica
        r_t = np.zeros(Nt, dtype=complex)

        # Somatório das N componentes multipercurso
        for n in range(self.N_multi):
            # Atraso do pulso para o percurso n
            t_deslocado = t - self.tau_n[n]
            s_atrasado = np.where(
                (t_deslocado >= atraso_ini)
                & (t_deslocado <= (atraso_ini + delta_t)),
                1.0,
                0.0,
            )

            # Fase variante no tempo
            fase_tempo = -1j * (
                self.phi_bar_n[n] - 2.0 * np.pi * self.nu_n[n] * t
            )

            r_t += alpha_n[n] * np.exp(fase_tempo) * s_atrasado

        return t, s_t, r_t

    # Método de plotagem adaptado para suportar subeixos (ax)
    def plot_sinal_recebido(self, delta_t, ax=None, atraso_ini=0.0, Nt=50000):
        show_plot = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(9, 5))
            show_plot = True

        t, s_t, r_t = self.simular_sinal_rx(delta_t, atraso_ini=atraso_ini, Nt=Nt)

        # Módulo do sinal transmitido e recebido
        ax.plot(
            t, np.abs(s_t), color="#1f77b4", linewidth=1.8, label="Transmissão (TX)"
        )
        ax.plot(
            t, np.abs(r_t), color="#d62728", linewidth=1.5, alpha=0.9, label="Recepção (RX)"
        )

        # Escala horizontal fixa em múltiplos de delta_t
        ticks_t = np.linspace(0.0, 5.0 * delta_t, 6)
        labels_t = ["0", r"$\delta t$", r"$2\delta t$", r"$3\delta t$", r"$4\delta t$", r"$5\delta t$"]
        
        ax.set_xticks(ticks_t)
        ax.set_xticklabels(labels_t)
        ax.set_xlim([0, 5.0 * delta_t])
        
        # Limite Y padronizado para comparação justa
        ax.set_ylim([0, 1.4])

        B_w = 1.0 / delta_t
        sigma_tau_ns = self.espa_atraso * 1e9
        
        # Título interno mais limpo
        titulo = rf"$\delta t = 10^{{{int(np.log10(delta_t))}}}\text{{ s}} \quad|\quad B_w = 10^{{{int(np.log10(B_w))}}}\text{{ Hz}}$"
        ax.set_title(titulo, fontsize=11)

        ax.grid(True, linestyle="--", alpha=0.5)
        
        # Minimizar poluição visual nos subplots matriciais
        if show_plot:
            ax.set_xlabel(r"Tempo absoluto — $t$ (s)", fontsize=12)
            ax.set_ylabel(r"$|\tilde{r}(t)|$", fontsize=12)
            ax.legend(loc="upper right")
            plt.tight_layout()
            plt.show()


# --- Teste reproduzindo os cenários dos slides 32 a 36 em Matriz ---
if __name__ == "__main__":
    
    # Valores de espalhamento de atraso em nanosegundos [10 ns, 100 ns, 1000 ns]
    valores_ds_ns = [10, 100, 1000]
    
    # Valores de duração do pulso em segundos [100 ns, 10 us, 1 ms]
    larguras_pulso_s = [1e-7, 1e-5, 1e-3]

    # Criação do painel 3x3
    fig, axes = plt.subplots(3, 3, figsize=(15, 10))
    fig.suptitle("Impacto Conjunto do Espalhamento de Atraso e Largura do Pulso", fontsize=16, fontweight='bold', y=0.98)

    print("\n--- Simulação Matricial do Sinal Recebido ---")

    for i, ds_ns in enumerate(valores_ds_ns):
        ds_segundos = ds_ns * 1e-9
        
        # Instancia um novo canal forçando o Espalhamento de Atraso específico
        canal_rx = P_RX(ambiente="UMi", freq=3.5, N_multi=100, override_espa_atraso=ds_segundos)
        
        print(f"Linha {i+1} | Espalhamento: {ds_ns} ns | Condição: {'LoS' if canal_rx.flag_Loss == 1 else 'NLoS'}")
        
        # Plota os três pulsos diferentes para o mesmo canal (mesma linha)
        for j, delta_t in enumerate(larguras_pulso_s):
            ax = axes[i, j]
            canal_rx.plot_sinal_recebido(delta_t=delta_t, ax=ax)
            
            # Adiciona rótulo de linha (Espalhamento) apenas na primeira coluna
            if j == 0:
                ax.set_ylabel(rf"$\mathbf{{\sigma_\tau = {ds_ns}\text{{ ns}}}}$" + "\n\n$|\\tilde{r}(t)|$", fontsize=12)
            
            # Adiciona rótulo de coluna (Tempo) e legenda apenas onde faz sentido visual
            if i == 2:
                ax.set_xlabel(r"Tempo ($t$)", fontsize=12)
            if i == 0 and j == 2:
                ax.legend(loc="upper right", framealpha=0.9, fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 0.96], h_pad=2.0, w_pad=1.5)
    plt.show()