import plotly
import numpy as np
import plotly
import plotly.graph_objs as go
N = 100
random_x = np.linspace(0, 1, N)
random_y2 = np.random.randn(N)

plotly.offline.plot({ "data": [go.Scatter(x=[1, 2, 3, 4], y=[4, 3, 2, 1])], "layout": go.Layout(title="hello world") })

#data = go.Scatter(  x = random_x,
#                    y = random_y2,
 #                    mode = 'lines',
 #                    name = 'lines')

#plotly.offline.plot(data, filename='scatter-mode.html')







