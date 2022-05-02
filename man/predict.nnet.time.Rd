\name{predict.nnet.time}
\alias{predict.nnet.time}

\title{
Prediction from a Survival Neural Network
}
\description{
Predict the survival of new observations based on a neural network estimated by using \code{nnet.time}.
}

\usage{
\method{predict}{nnet.time}(object, ..., newdata, newtimes)
}

\arguments{
  \item{object}{An object returned by the function \code{nnet.time}.}
  \item{...}{Further arguments passed.}
  \item{newdata}{An optional data frame containing covariate values at which to produce predicted values. There must be a column for every covariate included in \code{cov.quanti} and \code{cov.quali} included in the training sample. The default value is\code{NULL}, the predicted values are computed for the subjects of the training sample.}
  \item{newtimes}{The times at which to produce predicted values. The default value is\code{NULL}, the predicted values are computed for the observed times in the training data frame.}
}

\value{
\item{times}{A vector of numeric values with the times of the \code{predictions}.}
\item{predictions}{A matrix with the predictions of survivals of each subject (lines) for each observed times (columns).}
}


\author{
Yohann Foucher <Yohann.Foucher@univ-nantes.fr>

Camille Sabathe <camille.sabathe@univ-nantes.fr>
}

\examples{

# On a linux machine

#Install Python 

#sudo apt install python3-pip

#Install Python Miniconda for CPU computaion on a linux machine

#Download the file from https://docs.conda.io/en/latest/miniconda.html#linux-installers

#chmod +x Miniconda3-latest-Linux-x86_64.sh 
#./Miniconda3-latest-Linux-x86_64.sh 
#conda create --name torch python=3.8
#conda activate torch
#conda install pytorch torchvision torchaudio cpuonly -c pytorch
#conda install -c conda-forge pycox

#export PATH="/home/username/.local/bin:$PATH"

# remotes::install_github("RaphaelS1/survivalmodels")
# install_pycox(method = "auto",  conda = "auto",  pip = TRUE,  install_torch = TRUE )
# install_keras(method = "auto",  conda = "auto",  pip = TRUE,  install_tensorflow = TRUE)

# use_condaenv(condaenv = "~/miniconda3/envs/torch", conda = "auto", required = FALSE)

#data(dataDIVAT2)

# The estimation of the training model
#model<-nnet.time(times="times", failures="failures", data=dataDIVAT2,
#  cov.quanti=c("age"),  cov.quali=c("hla", "retransplant", "ecd"), n.nodes=10, decay=0,
#  batch.size=128, epochs=10)

# Predicted survival from the validation sample
#pred <- predict(model,
#  newdata=data.frame(age=c(52,52), hla=c(0,1), retransplant=c(1,1), ecd=c(0,1)))

#plot(y=pred$predictions[1,], x=pred$times, xlab="Time (years)", ylab="Predicted survival",
#     col=1, type="l", lty=1, lwd=2, ylim=c(0,1))

#lines(y=pred$predictions[2,], x=pred$times, col=2, type="l", lty=1, lwd=2)

#legend("bottomright", col=c(1,2), lty=1, lwd=2, c("Subject #1", "Subject #2"))
}

\keyword{Predicting}