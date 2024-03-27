import math
import numpy as np
from minisom import MiniSom


class Som:

    def __init__(self, features: list, labels: list, sigma=1.5, lr=0.5):
        self.labels = labels
        self.features = features

        self.sigma = sigma
        self.lr = lr
        self.random_seed = 123

        self.frequencies = None
        self.distance_map = None
        self.winning_neurons = None
        self.map_height = None
        self.map_width = None
        self.fn_get_label_marker = None
        self.fn_get_label_color = None

        self.som = None

    def auto_size(self):
        map_size = (10 * math.sqrt(len(self.labels)))
        map_height = map_width = math.ceil(math.sqrt(map_size))

        self.map_height = map_height
        self.map_width = map_width



        return self.map_height, self.map_width

    def run(self, num_iteration=1000):
        n_features = len(self.features[0])

        if self.map_width is None:
            self.auto_size()

        print(f'Number of features: {n_features}')
        print(f'(map_height, map_width) = ({self.map_height}, {self.map_width})')

        self.som = MiniSom(x=self.map_height, y=self.map_width, input_len=n_features,
                           sigma=self.sigma, learning_rate=self.lr,
                           neighborhood_function='gaussian',
                           random_seed=self.random_seed)

        self.som.pca_weights_init(self.features)
        self.som.train(data=self.features, num_iteration=num_iteration, verbose=True)  # random training

        self.frequencies = self.som.activation_response(self.features)
        self.distance_map = self.som.distance_map()
        print(f'Frequencies:\n {np.array(self.frequencies, np.uint)}')

        self.winning_neurons = np.array([self.som.winner(x) for x in self.features])

    def plot(self, offset_x=0.5, offset_y=0.5, marker_size=10):
        import matplotlib.pyplot as plt

        plt.figure(figsize=(self.map_height, self.map_width))

        u_matrix = self.distance_map.T
        plt.pcolor(u_matrix, cmap='bone_r')
        plt.colorbar()

        for i in range(len(self.features)):
            label = self.labels[i]
            w = self.winning_neurons[i]
            plt.plot(w[0] + offset_x, w[1] + offset_y,
                     self.fn_get_label_marker(label),
                     markeredgecolor=self.fn_get_label_color(label),
                     markerfacecolor='None', markersize=marker_size, markeredgewidth=1)

        return plt

    def plot_distance_map(self, ax, fig):
        """Plot the distance map"""
        p = ax.pcolor(self.distance_map.T, cmap='bone_r')  # cmap='Blues'
        # ax.colorbar()
        fig.colorbar(p, ax=ax)

    def plot_clusters_scatter(self, ax):
        """
        Create a scatter plot of the winning neurons.
        Each neuron is assigned the color of the cluster it belongs to.
        """

        # Add a random offset to avoid overlaps between points within the same cell
        offset = np.random.uniform(low=-0.4, high=0.4, size=(len(self.features), 2))
        winning_neurons = self.winning_neurons + offset

        # Define the colors based on the labels
        colors = ['#ff0400', 'g', '#e88325']
        label_colors = [self.fn_get_label_color(x) for x in self.labels]

        # Create the scatter plot
        # 1st column represent x and second, y coordinate
        ax.scatter(winning_neurons[:, 0], winning_neurons[:, 1], s=10, c=label_colors)

    def plot_clusters_markers(self, ax, marker_size=10):
        """
        Plot the winning neurons as markers.
        Each marker is assigned the color of the cluster ir belongs to.
        """
        for i, l in enumerate(self.labels):
            w = self.winning_neurons[i]
            marker = self.fn_get_label_marker(l)
            color = self.fn_get_label_color(l)
            ax.plot(w[0] + 0.5, w[1] + 0.5,
                    marker,
                    markeredgecolor=color,
                    markerfacecolor='None', markersize=marker_size, markeredgewidth=1)

        # legend
        """
        ax.legend(handles=[plt.Line2D([], [], color='#ff0400', marker='o', linestyle='None', label='Setosa'),
                           plt.Line2D([], [], color='green', marker='s', linestyle='None', label='Versicolor'),
                           plt.Line2D([], [], color='#e88325', marker='^', linestyle='None', label='Virginica')],
                  bbox_to_anchor=(1.5, 1.03))
        """

    def pp(self):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))

        self.plot_distance_map(axes[0], fig)
        self.plot_clusters_scatter(axes[1])
        self.plot_clusters_markers(axes[2])

        plt.suptitle("Plants species clusters")
        return plt
