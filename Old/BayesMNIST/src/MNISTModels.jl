#dimensions and parameters
const N=100
const D=784
const K=10
function modelmlp()
  imgs = MNIST.images()
  labels = MNIST.labels()

  #make each image into a feature vector
  #concatenate each feature vector into a matrix
  preX =hcat(float.(reshape.(imgs,:))...)
  X::CuMatrix{Float32} = preX |> gpu

  #creates a one-hot matrix of 0 vectors
  Y::CuMatrix{Float32} = onehotbatch(labels, 0:9)

  #dense layers
  #takes in all data for a matrix, creates an intermediate linear layer
  m = Chain(
    Dense(28^2, 32, relu),
    Dense(32, 10),
    softmax) |> gpu

  #some wacky loss function
  loss(x, y) = crossentropy(m(x), y)


  #repeat the sample by 200 times
  dataset = Base.Iterators.repeated((X, Y), 200)

  #callback to show the current loss
  evalcb = () -> @show(loss(X, Y))
  opt = ADAM()

  #train everything, calls the callback every 10 seconds (I think)
  @epochs 10 Flux.train!(loss, Flux.params(m), dataset, opt, cb = throttle(evalcb, 10))

  #evaluate the accuracy
  #checks how often two sets of data match
  accuracy(x, y) = mean(onecold(cpu(m(x))) .== onecold(cpu(y)))

  println("in sample accuracy")
  @show accuracy(X,Y)

  println("oos sample accuracy")
  tX = hcat(float.(reshape.(MNIST.images(:test), :))...) |> gpu
  tY = onehotbatch(MNIST.labels(:test), 0:9) |> gpu

  @show accuracy(tX, tY)

end
#=function modelbayes()
  imgs = MNIST.images()
  labels = MNIST.labels()

  #make each image into a feature vector
  #concatenate each feature vector into a matrix
  preX =hcat(float.(reshape.(imgs,:))...)
  X::CuMatrix{Float32} = preX |> gpu

  #creates a one-hot matrix of 0 vectors
  Y::CuMatrix{Float32} = onehotbatch(labels, 0:9)

  #NOTE: Now what? apply the VI step?
  #dense layers
  #takes in all data for a matrix, creates an intermediate linear layer
  m = Chain(
    Dense(D, K, hardtanh),
    Dense(32, 10),
    softmax) |> gpu

  @model bayesnetwork(x) = begin
    w ~ MvNormal(zeros(K), ones(K))
    b ~ MvNormal(zeros(K), ones(K))
    preds = m(x) + b
    y ~ Categorical(x*w + b)
  end

  Flux.Optimize.update()

  #some wacky loss function
  loss(x, y) = kldivergence(m(x), y)


  #repeat the sample by 200 times
  dataset = Base.Iterators.repeated((X, Y), 5000)

  #callback to show the current loss
  evalcb = () -> @show(loss(X, Y))
  opt = ADAM()

  #train everything, calls the callback every 10 seconds (I think)
  @epochs 1 Flux.train!(loss, Flux.params(m), dataset, opt, cb = throttle(evalcb, 10))

  #evaluate the accuracy
  #checks how often two sets of data match
  accuracy(x, y) = mean(onecold(cpu(m(x))) .== onecold(cpu(y)))

  println("in sample accuracy")
  @show accuracy(X,Y)

end=#
#=function modelcnn()
  imgs = MNIST.images()
  labels = MNIST.labels()

  #make each image into a feature vector
  #concatenate each feature vector into a matrix
  preX =hcat(float.(reshape.(imgs,:))...)
  X::CuMatrix{Float32} = preX |> gpu

  #creates a one-hot matrix of 0 vectors
  Y::CuMatrix{Float32} = onehotbatch(labels, 0:9)

  #make minibatches
  function makeminibatch(X, Y, idxs)
    Xbatch = Array{Float32}(undef, size(X[1])..., 1, length(idxs))
    for i in 1:length(idxs)
        Xbatch[:, :, :, i] = Float32.(X[idxs[i]])
    end
    Ybatch = onehotbatch(Y[idxs], 0:9)
    return (Xbatch, Ybatch)
  end

  batchindices = Base.iterators.partition(1:length(imgs),N)
  trainset = [makeminibatch(imgs, labels, i) for i in batchindices]

  m = Chain(
    Conv((3))


  @model bayesnetwork(x) = begin
    w ~ MvNormal(zeros(K), ones(K))
    b ~ MvNormal(zeros(K), ones(K))
    y ~ Categorical(x*w + b)
  end

  #some wacky loss function
  loss(x, y) = kldivergence(m(x), y)


  #repeat the sample by 200 times
  dataset = Base.Iterators.repeated((X, Y), 5000)

  #callback to show the current loss
  evalcb = () -> @show(loss(X, Y))
  opt = ADAM()

  #train everything, calls the callback every 10 seconds (I think)
  @epochs 10 Flux.train!(loss, Flux.params(m), dataset, opt, cb = throttle(evalcb, 10))

  #evaluate the accuracy
  #checks how often two sets of data match
  accuracy(x, y) = mean(onecold(cpu(m(x))) .== onecold(cpu(y)))

  println("in sample accuracy")
  @show accuracy(X,Y)

end=#
