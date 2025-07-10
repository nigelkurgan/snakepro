#This script is an adaption from alphaStats PCA class
############################################### IMPORTS ######################################
import pandas as pd
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

from sklearn.preprocessing import StandardScaler
from sklearn.covariance import MinCovDet
from sklearn.decomposition import PCA

from scipy.stats import chi2

######################################### PCA FUCTIONS ###################################

class RobustPCA:
    def __init__(self, dataset, group,  circle,  output_dir):
        self.dataset = dataset
        self.circle = circle
        self.group = group
        self.plot = None
        self.plot_distance = None
        self.plot_M_distance = None
        self.output_dir = output_dir

        sample_names, group_color = self._prepare_df()

        self._pca()

        self._plot(sample_names=sample_names, group_color=group_color)

        self._plot_distances(sample_names=sample_names, group_color=group_color)
        
        self._plot_M_distances(sample_names=sample_names, group_color=group_color)


    def _add_circles_to_scatterplot(self, fig):
        fig_dict = fig.to_plotly_json().get("data")
        for group in fig_dict:
            
            x_vector = group.get("x")
            y_vector = group.get("y")
            
            group_color = group.get("marker").get("color")
            fig.add_shape(type="circle",xref="x",yref="y",x0=min(x_vector),y0=min(y_vector),x1=max(x_vector),y1=max(y_vector),
                opacity=0.2,fillcolor=group_color,line_color=group_color,)
        return fig

    def _prepare_df(self):
        if self.group:
            mat = self.dataset._subset()
            self.dataset.metadata[self.group] = self.dataset.metadata[self.group].apply(str)
            group_color = self.dataset.metadata[self.group]
            sample_names = self.dataset.metadata[self.dataset.sample].to_list()
        else:
            mat = self.dataset.mat
            group_color = self.group
            sample_names = mat.reset_index(level=0)["index"].to_list()
        self.prepared_df = mat

        return sample_names, group_color

    def _pca(self):
        pca = PCA()
        self.components = pca.fit_transform(self.prepared_df)       
        self.labels = {"0": "PC 1 (%.2f%%)" % (pca.explained_variance_ratio_[0] * 100),
            "1": "PC 2 (%.2f%%)" % (pca.explained_variance_ratio_[1] * 100),}

    def _plot_distances(self, sample_names, group_color):
        """
        Compute pairwise Euclidean distances between PCA scores, and annotate
        only points that are 'yellow' to 'red' in the Jet scale (i.e., above ~75% of the data range).
        Report the top 5 samples with the highest Euclidean distance.
        """
        # Compute Euclidean distances using the first 5 PCs
        euclidean = np.zeros(self.prepared_df.shape[0])
        for i in range(5):
            euclidean += (self.components[:, i] - np.mean(self.components[:, :5]))**2 / np.var(self.components[:, :5])

        components = pd.DataFrame(self.components[:, :5],
                                columns=["PC1", "PC2", "PC3", "PC4", "PC5"])
        components[self.dataset.sample] = sample_names
        components["Euclidean Distance"] = euclidean

        # Determine the top 5 samples (for summary)
        top5 = components.sort_values("Euclidean Distance", ascending=False).head(5)
        top5_samples = top5[self.dataset.sample].tolist()
        top5_string = f"Samples {', '.join(top5_samples)} are the top 5 with highest Euclidean distance."
        print(top5_string)
        self.top5_euclidean = top5_string

        # Compute distance range
        dist_min = components["Euclidean Distance"].min()
        dist_max = components["Euclidean Distance"].max()
        fraction_for_yellow = 0.60  # adjust as needed (e.g. 0.80 or 0.85)
        threshold_euclid = dist_min + fraction_for_yellow * (dist_max - dist_min)

        # Annotate only if above threshold
        components["annot_text"] = components.apply(
            lambda row: row[self.dataset.sample] if row["Euclidean Distance"] > threshold_euclid else "",
            axis=1
        )

        # Create the Plotly scatter plot
        fig = px.scatter(
            components,
            x="PC1", y="PC2",
            labels={"PC1": self.labels["0"], "PC2": self.labels["1"], "Euclidean Distance": "Distance"},
            color="Euclidean Distance",
            template="simple_white+alphastats_colors",
            color_continuous_scale="Jet",
            text="annot_text"
        )
        fig.update_traces(textposition='top center')
        fig.update_layout(showlegend=False, title="PCA with Euclidean Distance")
        fig.update_coloraxes(showscale=False)
        fig.write_image(f"{self.output_dir}/pca_euclidean_distance.png")
        self.plot_distance = fig

    def _plot_M_distances(self, sample_names, group_color):
        """
        Compute pairwise Mahalanobis distances between PCA scores, and annotate
        only points that are 'yellow' to 'red' in the Jet scale (i.e., above ~75% of the data range).
        Report the top 5 samples with the highest Mahalanobis distance.
        """
        from sklearn.covariance import MinCovDet
        robust_cov = MinCovDet().fit(self.components[:, :5])
        mahalanobis_distances = robust_cov.mahalanobis(self.components[:, :5])

        components = pd.DataFrame(self.components[:, :5],
                                columns=["PC1", "PC2", "PC3", "PC4", "PC5"])
        components[self.dataset.sample] = sample_names
        components["Mahalanobis Distance"] = mahalanobis_distances

        # Determine the top 5 samples for summary
        top5 = components.sort_values("Mahalanobis Distance", ascending=False).head(5)
        top5_samples = top5[self.dataset.sample].tolist()
        top5_string = f"Samples {', '.join(top5_samples)} are the top 5 with highest Mahalanobis distance."
        print(top5_string)
        self.top5_mahalanobis = top5_string

        # Compute distance range
        dist_min = components["Mahalanobis Distance"].min()
        dist_max = components["Mahalanobis Distance"].max()
        fraction_for_yellow = 0.60  # adjust as needed (e.g. 0.80 or 0.85)
        threshold_maha = dist_min + fraction_for_yellow * (dist_max - dist_min)

        # Annotate only if above threshold
        components["annot_text"] = components.apply(
            lambda row: row[self.dataset.sample] if row["Mahalanobis Distance"] > threshold_maha else "",
            axis=1
        )

        # Create the Plotly scatter plot
        fig = px.scatter(
            components,
            x="PC1", y="PC2",
            labels={"PC1": self.labels["0"], "PC2": self.labels["1"], "Mahalanobis Distance": "Distance"},
            color="Mahalanobis Distance",
            hover_data=[components[self.dataset.sample]],
            template="simple_white+alphastats_colors",
            color_continuous_scale="Jet",
            text="annot_text"
        )
        fig.update_traces(textposition='top center')
        fig.update_coloraxes(showscale=False)
        fig.update_layout(showlegend=False, title="PCA with Mahalanobis Distance")
        fig.write_image(f"{self.output_dir}/pca_M_distance.png")
        self.plot_M_distance = fig


    def _plot(self,sample_names, group_color):
        components =  pd.DataFrame(self.components)
        components[self.dataset.sample] = sample_names

        fig = px.scatter(components,x=0,y=1,labels=self.labels,color=group_color,
            hover_data=[components[self.dataset.sample]],template="simple_white+alphastats_colors")

        fig_dict = fig.to_plotly_json()
        data = fig_dict.get("data")

        for count, d in enumerate(data):
            hover = d.get("hovertemplate").replace("hover_data_0", "sample")
            fig_dict["data"][count]["hovertemplate"] = hover
        fig = go.Figure(fig_dict)

        if self.circle is True and self.group is not None:
            fig = self._add_circles_to_scatterplot(fig)

        if self.group:
            fig.update_layout(legend_title_text=self.group)

        fig.write_image(f"{self.output_dir}/pca_plot.png")
        self.plot = fig

