// update 1.0
pipeline {
agent {label 'merlin'}
    stages {
        stage ("build-aws-gatk4") {
            steps {
                 dir("wspace-gatk4") {
                 checkout scm
//                    checkout([$class: 'GitSCM',
//                    branches: [[name: '*/release']],
//                    gitTool: 'Default', 
//                    userRemoteConfigs: [[url: 'git@github.com:falcon-computing/gatk4.git']],
//                    extensions: [[$class: 'CloneOption', timeout: 120]]
//                        ])
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
                                version= sh(returnStdout: true, script: 'git describe --tag').trim()
                                sh "echo $version"
                                sh "cd ~/falcon2/tools/package; mv GATK4.jar GATK4-$version-aws.jar"
                                sh "cd ~/falcon2/tools/package; echo s3://fcs-cicd-test/release/aws/gatk4/GATK4-$version-aws.jar > latest"
//                                link = sh(returnStdout: true, script: 'cd ~/falcon2/tools/package; link=s3://fcs-cicd-test/release/aws/gatk4/GATK4-$version-aws.jar; echo $link; echo $link > latest')
                        	    sh "cd ~/falcon2/tools/package; aws s3 cp GATK4-$version-aws.jar s3://fcs-cicd-test/release/aws/gatk4/GATK4-$version-aws.jar"
                        	    sh "cd ~/falcon2/tools/package; aws s3 cp latest s3://fcs-cicd-test/release/aws/gatk4/latest"
                        	    sh "cd ~/falcon2/tools/package; rm -f latest"
                            }
                        }
                    }
                }
            }
        }


  
