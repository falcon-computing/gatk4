// update 1.0
pipeline {
agent {label 'merlin'}
    stages {
        stage ("build-local-gatk4") {
            steps {
                 dir("ws-gatk4") {
                    checkout([$class: 'GitSCM',
                    branches: [[name: '*/release']],
                    gitTool: 'Default', 
                    userRemoteConfigs: [[url: 'git@github.com:falcon-computing/gatk4.git']],
                    extensions: [[$class: 'CloneOption', timeout: 120]]
                        ])
                     script {
                        dir("falcon"){
                            sh "./gradlew clean install -Prelease &>> ../build.log --no-daemon"
                            }
                                version = sh(returnStdout: true, script: 'cat build.log | grep "Version:" | awk \'{print $2}\'')
                                sh "./gradlew bundle -Dfalcon.version=$version &>> build.log"
                                sh "mkdir -p export"
                                sh "cp build/libs/gatk.jar ./export/"
                                sh "mv ./export/gatk.jar ./export/GATK4.jar"
                                sh "mv ./export/GATK4.jar /curr/limark/falcon2/tools/package/GATK4.jar"
                                sh "rm -f build.log"
                                link = sh(returnStdout: true, script: 'cd ~/falcon2/tools/package; link=s3://fcs-cicd-test/release/aws/gatk4/GATK4.jar; echo $link; echo $link > latest')
                        	    sh "cd ~/falcon2/tools/package; aws s3 cp GATK4.jar s3://fcs-cicd-test/release/aws/gatk4/GATK4.jar"
                        	    sh "cd ~/falcon2/tools/package; aws s3 cp latest s3://fcs-cicd-test/release/aws/gatk4/latest"
                        	    sh "cd ~/falcon2/tools/package; rm -f latest"
                            }
                        }
                    }
                }
            }
        post {
            always {

                emailext attachLog: true, body: "${currentBuild.currentResult}: Job ${env.JOB_NAME} build ${env.BUILD_NUMBER}\n More info at: ${env.BUILD_URL}console",
                    recipientProviders: [[$class: 'DevelopersRecipientProvider'], [$class: 'RequesterRecipientProvider']],
                    subject: "Jenkins Build ${currentBuild.currentResult}: Job ${env.JOB_NAME}",
                    to: 'roshantha@limarktech.com'

        }
    }
}

  
