

function runbnntest()
  # Hide sampling progress.
  Turing.turnprogress(false);

  # Use reverse_diff due to the number of parameters in neural networks.
  Turing.setadbackend(:reverse_diff)

  # Number of points to generate.
  N = 80
  M = round(Int, N / 4)
  Random.seed!(1234)

  # Generate artificial data.
  x1s = rand(M) * 4.5; x2s = rand(M) * 4.5;
  xt1s = Array([[x1s[i] + 0.5; x2s[i] + 0.5] for i = 1:M])
  x1s = rand(M) * 4.5; x2s = rand(M) * 4.5;
  append!(xt1s, Array([[x1s[i] - 5; x2s[i] - 5] for i = 1:M]))

  x1s = rand(M) * 4.5; x2s = rand(M) * 4.5;
  xt0s = Array([[x1s[i] + 0.5; x2s[i] - 5] for i = 1:M])
  x1s = rand(M) * 4.5; x2s = rand(M) * 4.5;
  append!(xt0s, Array([[x1s[i] - 5; x2s[i] + 0.5] for i = 1:M]))

  # Store all the data for later.
  xs = [xt1s; xt0s]
  ts = [ones(2*M); zeros(2*M)]

  # Plot data points.
  function plot_data()
      x1 = map(e -> e[1], xt1s)
      y1 = map(e -> e[2], xt1s)
      x2 = map(e -> e[1], xt0s)
      y2 = map(e -> e[2], xt0s)

      Plots.scatter(x1,y1, color="red", clim = (0,1))
      Plots.scatter!(x2, y2, color="blue", clim = (0,1))
  end

  plot_data()


  function unpack(nn_params::AbstractVector)
      W₁ = reshape(nn_params[1:6], 3, 2);
      b₁ = reshape(nn_params[7:9], 3)

      W₂ = reshape(nn_params[10:15], 2, 3);
      b₂ = reshape(nn_params[16:17], 2)

      Wₒ = reshape(nn_params[18:19], 1, 2);
      bₒ = reshape(nn_params[20:20], 1)
      return W₁, b₁, W₂, b₂, Wₒ, bₒ
  end

  function nn_forward(xs, nn_params::AbstractVector)
      W₁, b₁, W₂, b₂, Wₒ, bₒ = unpack(nn_params)
      nn = Chain(Dense(W₁, b₁, tanh),
                 Dense(W₂, b₂, tanh),
                 Dense(Wₒ, bₒ, σ))
      return nn(xs)
  end

  alpha = 0.09
  sig = sqrt(1.0 / alpha)

  # Specify the probabalistic model.
  @model bayes_nn(xs, ts) = begin
      # Create the weight and bias vector.
      nn_params ~ MvNormal(zeros(20), sig .* ones(20))

      # Calculate predictions for the inputs given the weights
      # and biases in theta.
      preds = nn_forward(xs, nn_params)

      # Observe each prediction.
      for i = 1:length(ts)
          ts[i] ~ Bernoulli(preds[i])
      end
  end

  N = 5000
  ch = sample(bayes_nn(hcat(xs...), ts), HMC(0.05, 4), N)
  theta = ch[:nn_params].value.data

  # Plot the data we have.
  plot_data()

  # Find the index that provided the highest log posterior in the chain.
  _, i = findmax(ch[:lp].value.data)

  # Extract the max row value from i.
  i = i.I[1]

  # Plot the posterior distribution with a contour plot.
  x_range = collect(range(-6,stop=6,length=25))
  y_range = collect(range(-6,stop=6,length=25))
  Z = [nn_forward([x, y], theta[i, :])[1] for x=x_range, y=y_range]
  contour!(x_range, y_range, Z)
end
